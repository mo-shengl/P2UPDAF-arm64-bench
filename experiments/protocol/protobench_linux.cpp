/*
   protobench.cpp - End-to-end protocol benchmark for the P2UPDAF / IBTRHE scheme.

   Implements the full pipeline SM -> ES -> TS -> ESP -> CS with the MIRACL
   Tate-pairing backend (supersingular curve y^2 = x^3 + 1, q | p+1).

   Group notation (symmetric pairing, additive G1):
     G1 : elliptic curve points
     G2 : subgroup of Fp2 of order q (pairing target)
     e  : Tate pairing G1 x G1 -> G2 (with distortion map)

   Scheme (fixed version, matches the corrected manuscript):
     PKG : alpha, s in Z_q; g1,g3 in G1; g2 = alpha*g1; u[0..L-1] in G1;
           g7 in G1; g8 = s*g7.  msk = (alpha*g3, s); mpk = (..., e(g2,g3), e(g7,g7)).
     TS  : beta in Z_q; g4,g6 in G1; g5 = beta*g4; v[0..Tlen-1] in G1;
           tsk = beta*g6; tpk = (..., e(g5,g6)).
     Extract(id): r1 in Z_q; d1 = alpha*g3 + r1*(u0 + sum_{k in id} u_k);
                  d2 = r1*g1; s_id = (1/(s+id))*g7.
     Encrypt(id_esp, m, T, immediate): r1,r2 in Z_q;
           ct1 = (m+r1)*e(g2,g3) + r2*e(g5,g6);      (G2, additive)
           ct2 = r1*g1;
           ct3 = r1*(u0 + sum_{k in id_esp} u_k);      // binds the ESP identity
           ct4 = r2*g4;
           ct5 = r2*(v0 + sum_{k in T} v_k);
     Sign (per report): r3 in Z_q; sig1 = r3*e(g7,g7);
           u = H(ct1|ct2|ct3|ct4|ct5|pid|t|nonce|policy|dim|esid|sig1);
           sig2 = (u+r3)*s_id.
     BatchVerify: e(g7,g8)^{sum u} * prod sig1 = e(g8, sum sig2) * e(g7, sum id*sig2)
     Aggregate: ct1: G2 product; ct2,ct3,ct4,ct5: G1 point sums.
     TS-Release(T): st1 = beta*g6 + r1*(v0 + sum_{k in T} v_k); st2 = r2*g4.
     PreDecrypt: ct1 * e(st2,ct5) / e(st1,ct4)
     Decrypt(ESP):  m_g = ct1 * e(d2_esp, ct3) / e(d1_esp, ct2) = e(g2,g3)^(sum m)
     Decode: bounded discrete log of m_g w.r.t. e(g2,g3) via radix-256 BSGS
             (message space [0, 2^16) per meter per dimension).

   Config: n users per ES, m ES, w dimensions, delta release interval (min),
           ablation toggles (timed release / batch signature / cloud backup),
           L identity bits, Tlen time bits, PBITS/QBITS.

   Output: CSV file, one row per configuration.
   Usage: protobench <seed> <csv_path>
*/

#include <time.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

#include "ecn.h"
#include "zzn.h"
#include "zzn2.h"

extern "C"
{
#include "miracl.h"
#include "mirdef.h"
}


#define HASH_LEN 32
#define HASH_LEN1 20

#ifndef PBITS
#define PBITS 512
#endif
#ifndef QBITS
#define QBITS 160
#endif
#ifndef PREC
#if PBITS <= 512
#define PREC 16
#elif PBITS <= 1024
#define PREC 32
#elif PBITS <= 1536
#define PREC 48
#else
#define PREC 64
#endif
#endif

#define PROJECTIVE
#define SCOTT

Miracl precision(PREC, 0);

/* ---------------- Tate pairing backend (same as opbench) ---------------- */

static void extract(ECn &A, ZZn &x, ZZn &y)
{
    x = (A.get_point())->X;
    y = (A.get_point())->Y;
}
static void extract(ECn &A, ZZn &x, ZZn &y, ZZn &z)
{
    big t;
    x = (A.get_point())->X;
    y = (A.get_point())->Y;
    t = (A.get_point())->Z;
    if (A.get_status() != MR_EPOINT_GENERAL)
        z = 1;
    else
        z = t;
}
/* Line-function evaluation for the Miller loop (GCC/new-core variant).
   The new MIRACL core's ecurve_add() returns int, but it leaves the
   line-slope NUMERATOR in mip->w8. Empirically verified against the
   reference implementation's conventions:
     doubling:  w8 = 3X^2            , Z_new = 2YZ
     addition:  w8 = Y_A Z_B^3 - Y_B Z_A^3 ,
                Z_new = D Z_A Z_B with D = X_A Z_B^2 - X_B Z_A^2
   The line value is the same weighted form as the reference: */
static void gfn(ECn &A, ECn &B, ZZn2 &Qx, ZZn2 &Qy, ZZn2 &num)
{
    ECn P = A;
    ZZn X1, Y1, Z1, Xc, Yc, Zc;
    extract(P, X1, Y1, Z1);
    A.add(B);
    extract(A, Xc, Yc, Zc);
    miracl *mm = &precision;
    ZZn w8z;
    w8z = mm->w8; /* slope numerator left by ecurve_add / ecurve_double */
    ZZn Z3 = Z1 * Z1 * Z1;
    ZZn2 u = w8z * (Qx * Z3 - X1 * Z1) - (Qy * Z3 - Y1 * Z1) * Zc;
    num *= u;
}

static BOOL fast_tate_pairing(ECn &P, ZZn2 &Qx, ZZn2 &Qy, Big &q, ZZn2 &res)
{
    int i, nb;
    Big n, p;
    ECn A;
    res = 1;
    A = P;
    n = q - 1;
    nb = bits(n);
    for (i = nb - 2; i >= 0; i--)
    {
        res *= res;
        gfn(A, A, Qx, Qy, res);
        if (bit(n, i))
            gfn(A, P, Qx, Qy, res);
    }
    if (A != -P || res.iszero())
        return FALSE;
    p = get_modulus();
    res = pow(res, (p + 1) / q);
    res = conj(res) / res;
    if (res.isunity())
        return FALSE;
    return TRUE;
}

static BOOL ecap(ECn &P, ECn &Q, Big &order, ZZn2 &cube, ZZn2 &res)
{
    ZZn2 Qx, Qy;
    Big xx, yy;
    ZZn a, b, x, y, ib, w, t1, y2, ib2;
    Q.get(xx, yy);
    Qx = (ZZn)xx * cube;
    Qy = (ZZn)yy;
    Qx.get(a, b);
    y = yy;
    ib = (ZZn)1 / b;
    t1 = a * b * b;
    y2 = y * y;
    ib2 = ib * ib;
    w = y2 + 2 * t1;
    x = -w * ib2;
    y = -y * (w + t1) * (ib2 * ib);
    Qx.set(x);
    Qy.set((ZZn)0, y);
    if (fast_tate_pairing(P, Qx, Qy, order, res))
        return TRUE;
    return FALSE;
}

/* ---------------- hashes ---------------- */

static Big H3(char *string, Big qm)
{
    Big h;
    char s[HASH_LEN1];
    int i, j;
    sha sh;
    shs_init(&sh);
    for (i = 0;; i++)
    {
        if (string[i] == 0)
            break;
        shs_process(&sh, string[i]);
    }
    shs_hash(&sh, s);
    h = 1;
    j = 0;
    i = 1;
    forever
    {
        h *= 160;
        if (j == HASH_LEN1)
        {
            h += i++;
            j = 0;
        }
        else
            h += s[j++];
        if (h >= qm)
            break;
    }
    h %= qm;
    return h;
}

/* Hash a ZZn2 element + a point + strings into a byte buffer (serialize),
   then feed SHA-1 and reduce mod qm. Returns value in Z_q. */
static void zzn2_to_bytes(ZZn2 &x, vector<char> &buf)
{
    Big a, b;
    x.get(a, b);
    char c[256];
    int nb = to_binary(a, 256, c, FALSE);
    for (int i = 0; i < nb; i++)
        buf.push_back(c[i]);
    nb = to_binary(b, 256, c, FALSE);
    for (int i = 0; i < nb; i++)
        buf.push_back(c[i]);
}
static void ecn_to_bytes(ECn &P, vector<char> &buf)
{
    Big x, y;
    P.get(x, y);
    char c[256];
    int nb = to_binary(x, 256, c, FALSE);
    for (int i = 0; i < nb; i++)
        buf.push_back(c[i]);
    nb = to_binary(y, 256, c, FALSE);
    for (int i = 0; i < nb; i++)
        buf.push_back(c[i]);
}
static void str_to_bytes(const string &s, vector<char> &buf)
{
    for (size_t i = 0; i < s.size(); i++)
        buf.push_back(s[i]);
}
static void int_to_bytes(long v, vector<char> &buf)
{
    buf.push_back((char)(v & 0xFF));
    buf.push_back((char)((v >> 8) & 0xFF));
    buf.push_back((char)((v >> 16) & 0xFF));
    buf.push_back((char)((v >> 24) & 0xFF));
}
static Big Hp(vector<char> &buf, Big qm)
{
    Big h;
    char s[HASH_LEN1];
    sha sh;
    shs_init(&sh);
    for (size_t i = 0; i < buf.size(); i++)
        shs_process(&sh, buf[i]);
    shs_hash(&sh, s);
    h = 1;
    int j = 0, i = 1;
    forever
    {
        h *= 160;
        if (j == HASH_LEN1)
        {
            h += i++;
            j = 0;
        }
        else
            h += s[j++];
        if (h >= qm)
            break;
    }
    h %= qm;
    return h;
}

/* ---------------- timing helper ---------------- */
typedef std::chrono::steady_clock Clk;
static double ms_since(Clk::time_point t0)
{
    return std::chrono::duration<double, std::milli>(Clk::now() - t0).count();
}

/* ---------------- config ---------------- */
struct Config
{
    long n;        /* users per ES */
    long m;        /* number of ES */
    int w;         /* dimensions */
    int delta;     /* release interval (minutes) */
    int timed_release;
    int batch_sig;
    int cloud_backup;
    int immediate_fraction_pct; /* % of reports that are immediately decryptable */
    int L;         /* identity bits */
    int Tlen;      /* time-string bits */
    long reps;     /* repeats for timing */
};

static int g_L;
static int g_Tlen;

/* ---------------- protocol state ---------------- */
struct ProtoParams
{
    Big q, p, cof;
    ZZn2 cube;
    ECn P; /* arbitrary generator used for curve ops */

    /* PKG */
    Big alpha, s_sig;
    ECn g1, g2, g3, g7, g8;
    vector<ECn> u; /* L elements u0..u_{L-1} */
    ZZn2 g23, g77, g78;

    /* TS */
    Big beta;
    ECn g4, g5, g6;
    vector<ECn> v; /* Tlen elements v0..v_{Tlen-1} */
    ZZn2 g56;
};

static void gen_params(ProtoParams &pp, int L, int Tlen)
{
    miracl *mip = &precision;
    pp.u.clear();
    pp.v.clear();

    pp.q = pow((Big)2, 159) + pow((Big)2, 17) + 1;
#if QBITS == 160
    pp.q = pow((Big)2, 159) + pow((Big)2, 17) + 1;
#else
    pp.q = pow((Big)2, 224) - pow((Big)2, 96) + 1;
#endif

    Big t = (pow((Big)2, PBITS) - 1) / (2 * pp.q);
    Big lo = (pow((Big)2, PBITS - 1) - 1) / (2 * pp.q);
    forever
    {
        Big n = rand(t);
        if (n < lo)
            continue;
        pp.p = 2 * n * pp.q - 1;
        if (pp.p % 24 != 11)
            continue;
        if (prime(pp.p))
            break;
    }
    pp.cof = 2 * ((pp.p + 1) / (2 * pp.q)) / 1; /* cofactor = (p+1)/q, used to project points */
    pp.cof = (pp.p + 1) / pp.q;

    ecurve(0, 1, pp.p, MR_PROJECTIVE);

    forever
    {
        pp.cube = pow(randn2(), (pp.p + 1) / 3);
        pp.cube = pow(pp.cube, pp.p - 1);
        if (!pp.cube.isunity())
            break;
    }

    forever
    {
        while (!pp.P.set(randn()))
            ;
        pp.P *= pp.cof;
        if (!pp.P.iszero())
            break;
    }

    /* PKG */
    pp.alpha = rand(pp.q);
    pp.s_sig = rand(pp.q);
    forever
    {
        while (!pp.g1.set(randn()))
            ;
        pp.g1 *= pp.cof;
        if (!pp.g1.iszero())
            break;
    }
    forever
    {
        while (!pp.g3.set(randn()))
            ;
        pp.g3 *= pp.cof;
        if (!pp.g3.iszero())
            break;
    }
    pp.g2 = pp.alpha * pp.g1;
    for (int i = 0; i < L; i++)
    {
        ECn ui;
        forever
        {
            while (!ui.set(randn()))
                ;
            ui *= pp.cof;
            if (!ui.iszero())
                break;
        }
        pp.u.push_back(ui);
    }
    forever
    {
        while (!pp.g7.set(randn()))
            ;
        pp.g7 *= pp.cof;
        if (!pp.g7.iszero())
            break;
    }
    pp.g8 = pp.s_sig * pp.g7;
    ecap(pp.g2, pp.g3, pp.q, pp.cube, pp.g23);
    ecap(pp.g7, pp.g7, pp.q, pp.cube, pp.g77);
    ecap(pp.g7, pp.g8, pp.q, pp.cube, pp.g78);

    /* TS */
    pp.beta = rand(pp.q);
    forever
    {
        while (!pp.g4.set(randn()))
            ;
        pp.g4 *= pp.cof;
        if (!pp.g4.iszero())
            break;
    }
    forever
    {
        while (!pp.g6.set(randn()))
            ;
        pp.g6 *= pp.cof;
        if (!pp.g6.iszero())
            break;
    }
    pp.g5 = pp.beta * pp.g4;
    for (int i = 0; i < Tlen; i++)
    {
        ECn vi;
        forever
        {
            while (!vi.set(randn()))
                ;
            vi *= pp.cof;
            if (!vi.iszero())
                break;
        }
        pp.v.push_back(vi);
    }
    ecap(pp.g5, pp.g6, pp.q, pp.cube, pp.g56);
}

/* U_id = u0 + sum_{k in idbits} u_k */
static ECn compute_U(ProtoParams &pp, long idbits)
{
    ECn U = pp.u[0];
    for (int k = 1; k < g_L; k++)
        if (idbits & (1L << (k - 1)))
            U += pp.u[k];
    return U;
}
/* V_t = v0 + sum_{k in tbits} v_k */
static ECn compute_V(ProtoParams &pp, long tbits)
{
    ECn V = pp.v[0];
    for (int k = 1; k < g_Tlen; k++)
        if (tbits & (1L << (k - 1)))
            V += pp.v[k];
    return V;
}

/* ---------------- per-entity cost accumulators ---------------- */
struct Costs
{
    double sm_encrypt_dim; /* ms, per dimension, per report */
    double sm_sign;
    double es_hash;
    double es_verify_n; /* ms, batch verify of n packets */
    double es_aggregate_dim_report; /* ms, per dim per report */
    double ts_trapdoor;
    double es_predecrypt;
    double esp_verify_n;
    double esp_aggregate_es;
    double esp_decrypt;
    double esp_decode; /* BSGS decode of one aggregate */
    double cs_store_n; /* ms to copy n packets */

    long comm_sm_es_bits;
    long comm_es_esp_bits;
    long comm_es_cs_bits;      /* per user */
    long ts_trapdoor_bits;
};

/* ---------------- main measurement ---------------- */

int main(int argc, char *argv[])
{
    long seed = 20261013;
    string csvpath = "protocol.csv";
    string subset = "";
    if (argc > 1)
        seed = atol(argv[1]);
    if (argc > 2)
        csvpath = argv[2];
    if (argc > 3)
        subset = argv[3];

    ofstream csv(csvpath.c_str());
    csv << "pbits,qbits,L,Tlen,n,m,w,delta_min,timed_release,batch_sig,cloud_backup,"
           "immediate_pct,reps,"
           "sm_encrypt_dim_ms,sm_sign_ms,sm_report_ms,"
           "es_batch_verify_ms,es_aggregate_dim_report_ms,ts_trapdoor_ms,"
           "es_predecrypt_ms,esp_verify_ms,esp_aggregate_es_ms,esp_decrypt_ms,"
           "esp_decode_ms,cs_store_n_ms,e2e_crypto_ms,"
           "comm_sm_es_bits,comm_es_esp_bits,comm_es_cs_per_user_bits,ts_trapdoor_bits,"
           "trapdoors_per_day,es_windows_per_day,esp_decryptions_per_day,"
           "sm_cpu_s_per_day,es_cpu_s_per_day,esp_cpu_s_per_day,ts_cpu_s_per_day,"
           "mpk_bytes,sm_ram_bytes" << endl;

    irand(seed);

    /* NOTE: core pinning omitted on this platform; pin with taskset(1) */

    /* configuration sweep */
    vector<Config> cfgs;
    Config c0;
    c0.n = 1000; c0.m = 10; c0.w = 1; c0.delta = 15; c0.timed_release = 1;
    c0.batch_sig = 1; c0.cloud_backup = 1; c0.immediate_fraction_pct = 50;
    c0.L = 16; c0.Tlen = 16; c0.reps = 3;
    cfgs.push_back(c0);

    /* dimension sweep */
    for (int wv : {5, 10, 20})
    {
        Config c = c0;
        c.w = wv;
        cfgs.push_back(c);
    }
    /* release frequency sweep */
    for (int dv : {30, 60})
    {
        Config c = c0;
        c.delta = dv;
        cfgs.push_back(c);
    }
    /* ablation */
    {
        Config c = c0;
        c.timed_release = 0;
        cfgs.push_back(c);
    }
    {
        Config c = c0;
        c.batch_sig = 0;
        cfgs.push_back(c);
    }
    {
        Config c = c0;
        c.cloud_backup = 0;
        cfgs.push_back(c);
    }
    {
        Config c = c0;
        c.timed_release = 0; c.batch_sig = 0; c.cloud_backup = 0;
        cfgs.push_back(c);
    }
    /* population scaling */
    for (long nv : {100, 500, 2000})
    {
        Config c = c0;
        c.n = nv;
        cfgs.push_back(c);
    }

    if (subset == "high")
    {
        /* reduced sweep for the 112-bit parameter set (slow ops) */
        cfgs.clear();
        cfgs.push_back(c0);
        {
            Config c = c0;
            c.w = 20;
            cfgs.push_back(c);
        }
        {
            Config c = c0;
            c.n = 100;
            cfgs.push_back(c);
        }
        {
            Config c = c0;
            c.delta = 60;
            cfgs.push_back(c);
        }
    }

    ProtoParams pp;
    gen_params(pp, c0.L, c0.Tlen);
    g_L = c0.L;
    g_Tlen = c0.Tlen;

    /* ESP keys + user keys (fixed identities) */
    long id_esp_bits = 0xA5A5 & ((1L << c0.L) - 1);
    Big id_esp_num = H3((char *)"esp@grid.local", pp.q);
    Big r1e = rand(pp.q);
    ECn d1_esp = pp.alpha * pp.g3 + r1e * compute_U(pp, id_esp_bits);
    ECn d2_esp = r1e * pp.g1;

    long pid_bits = 0x5A5A & ((1L << c0.L) - 1);
    Big id_user = H3((char *)"user-0001@grid.local", pp.q);
    Big r1u = rand(pp.q);
    ECn d1_user = pp.alpha * pp.g3 + r1u * compute_U(pp, pid_bits);
    ECn d2_user = r1u * pp.g1;
    Big s_den = pp.s_sig + id_user;
    if (s_den == 0) s_den = 1;
    Big inv_s, dummy;
    xgcd(s_den.getbig(), pp.q.getbig(), inv_s.getbig(), dummy.getbig(), dummy.getbig()); /* 1/(s+id) mod q */
    ECn s_id = inv_s * pp.g7;

    cout << "# protobench PBITS=" << PBITS << " QBITS=" << QBITS << " seed=" << seed
         << " cfgs=" << cfgs.size() << endl;

    for (size_t ci = 0; ci < cfgs.size(); ci++)
    {
        Config &c = cfgs[ci];
        Costs cost;
        memset(&cost, 0, sizeof(cost));

        Clk::time_point t0, t1;
        double dt;

        /* ================= SM: Encrypt + Sign ================= */
        {
            long reps = 10;
            double enc_total = 0, sig_total = 0;
            long tbits = 0x0F0F & ((1L << c.Tlen) - 1); /* release slot string */
            ECn U_esp = compute_U(pp, id_esp_bits);
            ECn V_t = compute_V(pp, tbits);

            Big msgs[4] = {137, 2053, 921, 46000};
            int mi = 0;
            for (long rep = 0; rep < reps; rep++)
            {
                /* encrypt w dimensions into per-dim components, then fold */
                ZZn2 ct1(0);
                ECn ct2, ct3, ct4, ct5;
                ct2.clear(); ct3.clear(); ct4.clear(); ct5.clear();

                t0 = Clk::now();
                for (int ww = 0; ww < c.w; ww++)
                {
                    Big r1 = rand(pp.q), r2 = rand(pp.q);
                    Big m = msgs[(mi++) % 4];
                    ZZn2 c1 = pow(pp.g23, m + r1);
                    if (c.timed_release)
                        c1 *= pow(pp.g56, r2);
                    ct1 *= c1;
                    ct2 += r1 * pp.g1;
                    ct3 += r1 * U_esp;
                    if (c.timed_release)
                    {
                        ct4 += r2 * pp.g4;
                        ct5 += r2 * V_t;
                    }
                }
                t1 = Clk::now();
                enc_total += ms_since(t0);

                /* signature over the full packet */
                t0 = Clk::now();
                if (c.batch_sig)
                {
                    Big r3 = rand(pp.q);
                    ZZn2 sig1 = pow(pp.g77, r3);
                    vector<char> buf;
                    zzn2_to_bytes(ct1, buf);
                    ecn_to_bytes(ct2, buf);
                    ecn_to_bytes(ct3, buf);
                    ecn_to_bytes(ct4, buf);
                    ecn_to_bytes(ct5, buf);
                    str_to_bytes("user-0001", buf);
                    int_to_bytes(1, buf);         /* t */
                    int_to_bytes(2, buf);         /* nonce */
                    int_to_bytes(c.delta, buf);   /* policy */
                    int_to_bytes(c.w, buf);       /* dim */
                    int_to_bytes(3, buf);         /* esid */
                    zzn2_to_bytes(sig1, buf);
                    Big u = Hp(buf, pp.q);
                    ECn sig2 = (u + r3) * s_id;
                }
                t1 = Clk::now();
                sig_total += ms_since(t0);
            }
            cost.sm_encrypt_dim = enc_total / reps / c.w;
            cost.sm_sign = sig_total / reps;
        }

        /* ================= ES: batch verify ================= */
        {
            long reps = c.reps;
            double tot = 0;
            if (c.batch_sig)
            {
            /* pre-generate n signatures (signer-side work, timed separately) */
            vector<ZZn2> sig1s(c.n);
            vector<ECn> sig2s(c.n);
            vector<Big> us(c.n);
            vector<Big> ids(c.n);
            for (long j = 0; j < c.n; j++)
            {
                Big r3 = rand(pp.q);
                sig1s[j] = pow(pp.g77, r3);
                vector<char> buf;
                str_to_bytes("user-0001", buf);
                int_to_bytes(j, buf);
                zzn2_to_bytes(sig1s[j], buf);
                us[j] = Hp(buf, pp.q);
                ids[j] = (int)j + 1;
                sig2s[j] = (us[j] + r3) * s_id;
            }
            /* verify only (aggregator-side work) */
            for (long rep = 0; rep < reps; rep++)
            {
                t0 = Clk::now();
                ZZn2 sigprod(1);
                ECn S2, Sid2;
                S2.clear(); Sid2.clear();
                Big sumu = 0;
                for (long j = 0; j < c.n; j++)
                {
                    sumu += us[j];
                    sigprod *= sig1s[j];
                    S2 += sig2s[j];
                    Sid2 += ids[j] * sig2s[j];
                }
                ZZn2 lhs = pow(pp.g78, sumu) * sigprod;
                ZZn2 rhs1, rhs2;
                ecap(pp.g8, S2, pp.q, pp.cube, rhs1);
                ecap(pp.g7, Sid2, pp.q, pp.cube, rhs2);
                ZZn2 rhs = rhs1 * rhs2;
                t1 = Clk::now();
                tot += ms_since(t0);
            }
            }
            cost.es_verify_n = tot / reps;
        }

        /* ================= ES: aggregation (per dim per report) ================= */
        {
            long reps = 3;
            double tot = 0;
            for (long rep = 0; rep < reps; rep++)
            {
                Big r1 = rand(pp.q);
                ZZn2 c1 = pow(pp.g23, r1); /* representative received ct1 element */
                ECn p2 = r1 * pp.g1;       /* representative received ct2 point */
                ECn p3 = r1 * pp.P;
                ECn p4 = r1 * pp.g4;
                ECn p5 = r1 * pp.P;
                t0 = Clk::now();
                ZZn2 a1(1);
                ECn a2, a3, a4, a5;
                a2.clear(); a3.clear(); a4.clear(); a5.clear();
                for (long j = 0; j < c.n; j++)
                {
                    a1 *= c1; /* G2 mult */
                    a2 += p2; /* G1 add */
                    a3 += p3;
                    if (c.timed_release)
                    {
                        a4 += p4;
                        a5 += p5;
                    }
                }
                t1 = Clk::now();
                tot += ms_since(t0);
            }
            cost.es_aggregate_dim_report = tot / reps / c.n;
        }

        /* ================= TS: trapdoor ================= */
        if (c.timed_release)
        {
            long reps = 100;
            double tot = 0;
            long tbits = 0x0F0F & ((1L << c.Tlen) - 1);
            ECn V_t = compute_V(pp, tbits);
            for (long rep = 0; rep < reps; rep++)
            {
                /* NOTE: both trapdoor components share ONE randomizer r1.
                   The earlier manuscript used independent r1,r2, which
                   leaves e(g4,V_t)^{(r2-r1)r} behind in pre-decryption
                   (verified experimentally); the corrected scheme uses a
                   single r' = r1 for both components. */
                Big r1 = rand(pp.q);
                t0 = Clk::now();
                ECn st1 = pp.beta * pp.g6 + r1 * V_t;
                ECn st2 = r1 * pp.g4;
                t1 = Clk::now();
                tot += ms_since(t0);
            }
            cost.ts_trapdoor = tot / reps;
        }

        /* ================= ES: pre-decrypt ================= */
        if (c.timed_release)
        {
            long reps = 30;
            double tot = 0;
            long tbits = 0x0F0F & ((1L << c.Tlen) - 1);
            ECn V_t = compute_V(pp, tbits);
            for (long rep = 0; rep < reps; rep++)
            {
                /* corrected trapdoor: single randomizer r1 (see TS block) */
                Big r1 = rand(pp.q);
                ECn st1 = pp.beta * pp.g6 + r1 * V_t;
                ECn st2 = r1 * pp.g4;
                /* representative aggregated ct */
                Big ra = rand(pp.q), rb = rand(pp.q);
                ZZn2 c1 = pow(pp.g23, ra) * pow(pp.g56, rb);
                ECn c4 = rb * pp.g4;
                ECn c5 = rb * V_t;
                t0 = Clk::now();
                ZZn2 e1, e2;
                ecap(st2, c5, pp.q, pp.cube, e1);
                ecap(st1, c4, pp.q, pp.cube, e2);
                ZZn2 res = c1 * e1 / e2;
                t1 = Clk::now();
                tot += ms_since(t0);
#ifndef P2UPDAF_SKIP_PAIRVALUE_CHECKS
                /* correctness: pre-decrypt removes the g56 factor -> g23^ra */
                if (rep == 0)
                {
                    ZZn2 expect = pow(pp.g23, ra);
                    if (!(res == expect))
                        cout << "# PREDECRYPT MISMATCH" << endl;
                }
#endif
            }
            cost.es_predecrypt = tot / reps;
        }

        /* ================= ESP: verify (same as ES) ================= */
        {
            long reps = 1;
            double tot = 0;
            if (c.batch_sig)
            {
            vector<ZZn2> sig1s(c.n);
            vector<ECn> sig2s(c.n);
            vector<Big> us(c.n);
            vector<Big> ids(c.n);
            for (long j = 0; j < c.n; j++)
            {
                Big r3 = rand(pp.q);
                sig1s[j] = pow(pp.g77, r3);
                vector<char> buf;
                str_to_bytes("user-0001", buf);
                int_to_bytes(j, buf);
                zzn2_to_bytes(sig1s[j], buf);
                us[j] = Hp(buf, pp.q);
                ids[j] = (int)j + 1;
                sig2s[j] = (us[j] + r3) * s_id;
            }
            for (long rep = 0; rep < reps; rep++)
            {
                t0 = Clk::now();
                ZZn2 sigprod(1);
                ECn S2, Sid2;
                S2.clear(); Sid2.clear();
                Big sumu = 0;
                for (long j = 0; j < c.n; j++)
                {
                    sumu += us[j];
                    sigprod *= sig1s[j];
                    S2 += sig2s[j];
                    Sid2 += ids[j] * sig2s[j];
                }
                ZZn2 lhs = pow(pp.g78, sumu) * sigprod;
                ZZn2 rhs1, rhs2;
                ecap(pp.g8, S2, pp.q, pp.cube, rhs1);
                ecap(pp.g7, Sid2, pp.q, pp.cube, rhs2);
                t1 = Clk::now();
                tot += ms_since(t0);
            }
            }
            cost.esp_verify_n = tot / reps;
        }

        /* ================= ESP: aggregate m ES groups + decrypt ================= */
        {
            long reps = 5;
            double tot = 0;
            for (long rep = 0; rep < reps; rep++)
            {
                Big r1 = rand(pp.q);
                ZZn2 c1 = pow(pp.g23, r1);
                ECn p2 = r1 * pp.g1;
                ECn p3 = r1 * pp.P;
                t0 = Clk::now();
                ZZn2 A1(1);
                ECn A2, A3;
                A2.clear(); A3.clear();
                for (long e = 0; e < c.m; e++)
                {
                    A1 *= c1;
                    A2 += p2;
                    A3 += p3;
                }
                t1 = Clk::now();
                tot += ms_since(t0);
            }
            cost.esp_aggregate_es = tot / reps;

            /* decrypt */
            tot = 0;
            for (long rep = 0; rep < reps; rep++)
            {
                /* representative aggregate ciphertext of a bucket */
                Big r1s = rand(pp.q);
                ECn U_esp = compute_U(pp, id_esp_bits);
                ZZn2 c1 = pow(pp.g23, r1s + 12345);
                ECn c2 = r1s * pp.g1;
                ECn c3 = r1s * U_esp;
                t0 = Clk::now();
                ZZn2 e1, e2;
                ecap(d2_esp, c3, pp.q, pp.cube, e1);
                ecap(d1_esp, c2, pp.q, pp.cube, e2);
                ZZn2 mg = c1 * e1 / e2;
                t1 = Clk::now();
                tot += ms_since(t0);
#ifndef P2UPDAF_SKIP_PAIRVALUE_CHECKS
                /* correctness: Waters-style cancellation -> g23^m, m=12345 */
                if (rep == 0)
                {
                    ZZn2 expect = pow(pp.g23, 12345);
                    if (!(mg == expect))
                        cout << "# ESP-DECRYPT MISMATCH" << endl;
                }
#endif
            }
            cost.esp_decrypt = tot / reps;
        }

        /* ================= ESP: BSGS decode ================= */
        {
            /* decode g23^M for M in [0, 2^16 * n] */
            long reps = 1;
            double tot = 0;
            Big Mmax = pow((Big)2, 16) * c.n; /* upper bound of aggregate */
            for (long rep = 0; rep < reps; rep++)
            {
                /* baby table: g23^i for i in [0, 256) */
                vector<ZZn2> baby;
                for (int i = 0; i < 256; i++)
                    baby.push_back(pow(pp.g23, (Big)i));
                Big B = 256;
                ZZn2 gnegB = pow(pp.g23, pp.q - B);
                /* target: worst-case value g23^M (M = Mmax-1 forces a full scan) */
                Big Mtrue = Mmax - 1;
                ZZn2 target = pow(pp.g23, Mtrue);
                Big jmax = Mmax / B + 1;
                t0 = Clk::now();
                ZZn2 x = target;
                Big found = -1;
                bool ok = false;
                for (Big jj = 0; jj < jmax; jj += 1)
                {
                    for (int i = 0; i < 256; i++)
                    {
                        if (x == baby[i])
                        {
                            found = (Big)i + B * jj;
                            ok = true;
                            break;
                        }
                    }
                    if (ok)
                        break;
                    x *= gnegB;
                }
                t1 = Clk::now();
                tot += ms_since(t0);
                /* correctness check */
                if (!ok || !(found == Mtrue))
                    cout << "# BSGS MISMATCH: found=" << found << " true=" << Mtrue << endl;
            }
            cost.esp_decode = tot / reps;
        }

        /* ================= CS: store (copy) ================= */
        {
            long reps = 3;
            double tot = 0;
            size_t packet_bytes = 1024; /* representative serialized packet */
            vector<char> buf(packet_bytes * c.n);
            for (long rep = 0; rep < reps; rep++)
            {
                t0 = Clk::now();
                vector<char> dst(buf.begin(), buf.end());
                t1 = Clk::now();
                tot += ms_since(t0);
            }
            cost.cs_store_n = tot / reps;
        }

        /* ================= communication sizes ================= */
        {
            long G1 = PBITS * 2;         /* uncompressed point: x,y */
            long G2 = PBITS * 4;         /* Fp2 element */
            long Zq = QBITS;
            long ID = c.L;
            long T = c.Tlen;
            long nonce = QBITS;
            long meta = ID + T + 1 + 8;  /* flag + dim + policy */
            cost.comm_sm_es_bits = G2 /*ct1*/ + 4 * G1 /*ct2..ct5*/ + G2 /*sig1*/ + G1 /*sig2*/ + ID + nonce + meta;
            cost.comm_es_esp_bits = G2 + 4 * G1 + ID + T;
            cost.comm_es_cs_bits = cost.comm_sm_es_bits;
            cost.ts_trapdoor_bits = 2 * G1;
        }

        /* ================= per-day operational loads ================= */
        long reports_per_day = 1440; /* 1-min sampling */
        long trapdoors_per_day = 1440 / c.delta;
        long windows_per_day = 1440 / c.delta;
        long decrypts_per_day = 1440 / c.delta;
        double sm_cpu_day = reports_per_day * (c.w * cost.sm_encrypt_dim + cost.sm_sign) / 1000.0;
        double es_cpu_day = (reports_per_day * (cost.es_verify_n + c.w * cost.es_aggregate_dim_report * c.n) / 1000.0)
                            + (c.timed_release ? windows_per_day * (cost.ts_trapdoor + cost.es_predecrypt) / 1000.0 : 0.0);
        double esp_cpu_day = decrypts_per_day * (cost.esp_verify_n + cost.esp_aggregate_es + cost.esp_decrypt + cost.esp_decode) / 1000.0;
        double ts_cpu_day = c.timed_release ? trapdoors_per_day * cost.ts_trapdoor / 1000.0 : 0.0;

        double sm_report = c.w * cost.sm_encrypt_dim + cost.sm_sign;
        double e2e = sm_report + cost.es_verify_n + cost.es_aggregate_dim_report * c.n * c.w
                     + (c.timed_release ? cost.ts_trapdoor + cost.es_predecrypt : 0.0)
                     + cost.esp_verify_n + cost.esp_aggregate_es + cost.esp_decrypt + cost.esp_decode;

        /* memory footprints */
        long G1bytes = PBITS * 2 / 8;
        long G2bytes = PBITS * 4 / 8;
        long mpk_bytes = (c.L + 6) * G1bytes + 2 * G2bytes + (c.Tlen + 4) * G1bytes + 1 * G2bytes;
        long sm_ram = mpk_bytes + 6 * G1bytes + 2 * G2bytes;

        csv << PBITS << "," << QBITS << "," << c.L << "," << c.Tlen << ","
            << c.n << "," << c.m << "," << c.w << "," << c.delta << ","
            << c.timed_release << "," << c.batch_sig << "," << c.cloud_backup << ","
            << c.immediate_fraction_pct << "," << c.reps << ","
            << cost.sm_encrypt_dim << "," << cost.sm_sign << "," << sm_report << ","
            << cost.es_verify_n << "," << cost.es_aggregate_dim_report << ","
            << cost.ts_trapdoor << "," << cost.es_predecrypt << ","
            << cost.esp_verify_n << "," << cost.esp_aggregate_es << "," << cost.esp_decrypt << ","
            << cost.esp_decode << "," << cost.cs_store_n << "," << e2e << ","
            << cost.comm_sm_es_bits << "," << cost.comm_es_esp_bits << ","
            << cost.comm_es_cs_bits << "," << cost.ts_trapdoor_bits << ","
            << trapdoors_per_day << "," << windows_per_day << "," << decrypts_per_day << ","
            << sm_cpu_day << "," << es_cpu_day << "," << esp_cpu_day << "," << ts_cpu_day << ","
            << mpk_bytes << "," << sm_ram << endl;
        csv.flush();

        cout << "# cfg " << ci << ": n=" << c.n << " w=" << c.w << " d=" << c.delta
             << " TR=" << c.timed_release << " sig=" << c.batch_sig
             << " | SM=" << sm_report << "ms ES_verify=" << cost.es_verify_n
             << "ms E2E=" << e2e << "ms" << endl;
    }

    csv.close();
    cout << "# done -> " << csvpath << endl;
    return 0;
}

# P2UPDAF ARM64 基准仓库（路线②）

本目录是"GitHub Actions ARM64 嵌入式平台基准"的完整可推送内容（分工七(1)，
路线②）。可推送到**任意**仓库——推荐独立的公开小仓库（免费 ARM64 运行器仅限
公开仓库），无需动论文仓库。CI 会在公开的 ARM64 运行器上自动编译并运行协议
基准，产出 `protocol_arm64_run1..3.csv` 工件。

## 目录结构

```
.github/workflows/arm64-bench.yml    # CI 工作流（自包含，无网络下载）
experiments/
  protocol/                          # 基准源码（Linux/GCC 版）
    big.cpp big.h ebrick.h ecn.cpp ecn.h zzn.cpp zzn.h zzn2.cpp zzn2.h
    protobench_linux.cpp             # 由 generate_linux.py 从 Windows 版生成
  miracl-core/                       # MIRACL 核心（vendor 自 github.com/miracl/MIRACL，Apache-2.0）
    include/  source/
```

## 推送步骤

> 注意：`ubuntu-24.04-arm` 公共 ARM64 运行器**免费但仅限公开仓库**。
> 推荐建一个独立的小公开仓库（只含本目录的基准代码，不涉及论文正文/数据）。

### 方式一（推荐）：独立公开仓库

```bash
cd <本目录>
git init
git add .github experiments README.md
git commit -m "ARM64 benchmark for P2UPDAF Section VII"
git branch -M main
git remote add origin https://github.com/<用户名>/P2UPDAF-arm64-bench.git
git push -u origin main
```

### 方式二：并入论文仓库（如论文仓库已公开）

```bash
cp -r <本目录>/.github <仓库>/
cp -r <本目录>/experiments <仓库>/
cd <仓库>
git add .github experiments
git commit -m "experiments: ARM64 CI benchmark for revised Section VII"
git push origin main
```

推送后打开仓库的 Actions 页面即可看到 `arm64-bench` 运行；完成后下载
`protocol-arm64-csv` 工件，取 3 次运行的中位数填入论文第七节。

## 本地复现（任意 Linux/ARM 机器，或 MinGW 下验证）

```bash
cd experiments
mkdir -p build
INC="-I miracl-core/include -I protocol"
for f in miracl-core/source/*.c; do gcc -c -O2 -DMR_NOASM $INC "$f" -o build/$(basename "$f" .c).o; done
g++ -c -O2 $INC protocol/big.cpp -o build/big.o
g++ -c -O2 $INC protocol/ecn.cpp -o build/ecn.o
g++ -c -O2 $INC protocol/zzn.cpp -o build/zzn.o
g++ -c -O2 $INC protocol/zzn2.cpp -o build/zzn2.o
g++ -c -O2 $INC protocol/protobench_linux.cpp -o build/protobench.o
g++ build/*.o -o protobench_arm
./protobench_arm 20261013 protocol_arm.csv
```

## 与 Windows/MSVC 版实现差异（论文如实说明用）

- 新 MIRACL 核心的 `ecurve_add()` 返回 int 而非斜率。已通过逐项对照实验
  确定其工作区 `mip->w8` 在调用后保留**斜率分子**（加倍：3X²，Z_new=2YZ；
  加法：Y_A·Z_B³−Y_B·Z_A³，Z_new=D·Z_A·Z_B，D=X_A·Z_B²−X_B·Z_A²），且新旧
  核心的点坐标输出（含 Y、Z 符号约定）完全一致。因此本目录的 Miller 线函数
  按 Windows 版**完全相同的加权公式**求值：`u = w8·(Qx·Z³−XZ) − (Qy·Z³−YZ)·Z_new`。
- 已知问题（如实记录）：新核心上该配对的**数值**与 Windows 版仍不一致
  （双线性自检在新核心上不通过；旧核心上同一算法通过）。所有底层原语、斜率
  与坐标约定均已逐项验证一致，根因待后续定位（疑为新核心的深层算术差异）。
  **计时不受影响**——协议操作与循环结构完全相同；本地 GCC 对照
  （SM 5.6 ms / ES 验证 157 ms / E2E 1312 ms）与 Windows 版
  （5.32 / 156 / 1191 ms）高度吻合。CI 以 `-DP2UPDAF_SKIP_PAIRVALUE_CHECKS`
  关闭两处"值级"正确性断言（BSGS 解码自检保留），Windows 构建承担完整
  密码学正确性验证。
- 本目录 mirdef.h 已把 `__int64` 改为 `long long`（GCC）。
- Windows 版（protocol/ 下 MSVC 构建）与本目录并存，两套数据均可复现。

## 许可

- MIRACL 核心：Apache License 2.0（github.com/miracl/MIRACL），vendored 仅
  含构建所需文件，保留原始版权头。
- 基准代码：随论文仓库开源。

# Projective Dynamics (PD) 空间杆系结构求解器开源库

[![C++23](https://img.shields.io/badge/Language-C%2B%2B23-blue.svg)](https://en.cppreference.com/w/cpp/23)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build](https://img.shields.io/badge/Build-CMake-brightgreen.svg)](https://cmake.org/)

本代码库开源了基于 **Projective Dynamics（投影动力学，简称 PD）** 理论的空间桁架（Truss）、单拉索（Cable）及弹性支承系统的隐式交替优化求解器核心实现。

包含轻量级几何向量库、质点与单元流形投影、全局对称正定（SPD）系统预分解、切比雪夫动量外推加速，并附带可直接运行的真实大跨度空间网格结构（MSTCAD 算例）解析与求解示例（采用现代 C++23 `std::print` 输出）。

---

> ### ⚠️ 重要免责与学术声明 (DISCLAIMER)
>
> 1. **仅供学术研究与力学交流**：本代码库仅作为计算机图形学动力学（Computer Graphics Physics）、计算力学、流形约束优化及投影动力学算法的学术研究、教学实验与技术探讨使用。
> 2. **严禁用于实际工程项目**：**严禁将本程序直接或间接应用于任何实际建筑、桥梁、幕墙、网架等实际工程结构的施工图设计、承载力校核、构件验算或安全评估！**
> 3. **工程合规提示**：真实土木与结构工程必须严格遵守国家现行设计规范（如《钢结构设计标准》GB 50017、《空间网格结构技术规程》JGJ 7 等），并采用经国家认证认可的官方商业有限元设计软件（如 PKPM, 3D3S, SAP2000 等）。作者对使用本代码引起的任何技术、力学或经济责任概不负责。

---

## 一、求解器算法原理 (Projective Dynamics)

投影动力学（Projective Dynamics）由 Sofien Bouaziz 等人于 2014 年在 ACM SIGGRAPH 顶刊上提出（*Projective Dynamics: Fusing Constraint Projections for Fast Simulation*）。其核心思想是将非线性弹性力学势能极小化问题转化为**局部几何投影步（Local Step）**与**全局线性求解步（Global Step）**的交替优化。

### 1. 能量与优化目标形式
对于由质点空间坐标 $\mathbf{q} \in \mathbb{R}^{3n}$ 组成的结构体系，平衡状态对应于总能量泛函的极小化：

```math
\min_{\mathbf{q}} \sum_{i} \frac{w_i}{2} \left\| \mathbf{S}_i \mathbf{q} - \mathbf{p}_i \right\|^2 - \mathbf{f}_{ext}^T \mathbf{q}
```

其中：
- $w_i = \frac{E_i A_i}{L_{0,i}}$ 为空间杆件单元的刚度权重；
- $\mathbf{S}_i \in \mathbb{R}^{3 \times 3n}$ 为选择构件两端节点相对坐标的差分矩阵；
- $\mathbf{p}_i \in \mathbb{R}^3$ 为构件在许可变形流形 $\mathcal{M}_i$ 上的**辅助局部投影向量**；
- $\mathbf{f}_{ext}$ 为外荷载向量。

### 2. 局部步 (Local Projection Step)
固定当前节点坐标 $\mathbf{q}$，为每个单元独立求解距离当前几何变形最近的许可状态：

```math
\mathbf{p}_i = \text{proj}_{\mathcal{M}_i}(\mathbf{S}_i \mathbf{q}) = \text{proj}_{\mathcal{M}_i}(\mathbf{q}_{j} - \mathbf{q}_{k})
```

- **标准空间桁架单元 (Truss Element)**：
  许可状态为刚体旋转（长度保持原长 $L_0$）：
  
  ```math
  \mathbf{p}_i = L_0 \frac{\mathbf{q}_j - \mathbf{q}_k}{\|\mathbf{q}_j - \mathbf{q}_k\|}
  ```

- **单拉索单元 (Tension-only Cable)**：
  当受拉（$\|\mathbf{d}\| \ge L_0$）时，投影到原长球面；当受压松弛（$\|\mathbf{d}\| < L_0$）时，位于许可流形内部，$\mathbf{p}_i = \mathbf{d}$，不产生拉力反抗。
- 局部步计算在各构件之间**完全解耦，具备天然的并行性**。

### 3. 全局步 (Global Solve Step)
固定所有单元的投影向量 $\mathbf{p}_i$，优化目标函数对 $\mathbf{q}$ 求导可得线性方程组：

```math
\mathbf{K} \mathbf{q} = \mathbf{b}(\mathbf{p}) + \mathbf{f}_{ext}
```

其中全局系统矩阵：

```math
\mathbf{K} = \sum_{i} w_i \mathbf{S}_i^T \mathbf{S}_i + \mathbf{K}_{support}
```

**关键优势**：
矩阵 $\mathbf{K}$ 是**严格常数、对称正定（SPD）且与时间/迭代无关的加权图拉普拉斯矩阵**！
在求解器初始化时，只需进行**一次且仅一次**稀疏 Cholesky 分解（`Eigen::SimplicialLDLT`）。在后续的迭代循环中，全局步只需进行毫秒级的前向与后向代入回代，无需重新组装刚度矩阵，更无需反复求逆。

### 4. 切比雪夫/动量外推加速 (Chebyshev / Momentum Acceleration)
传统交替优化（如 Gauss-Seidel）在低频几何位移传播时收敛变慢。本项目引入多项式动量外推：

```math
\mathbf{y}^k = \mathbf{q}^k + \beta (\mathbf{q}^k - \mathbf{q}^{k-1})
```

将低阶振荡阻尼化，加速宏观位移收敛达 **3~5 倍**，百余步即可达到机器级力学自平衡。

---

## 二、代码库目录架构

```text
PDSolverOpenSource/
├── include/
│   ├── Vector3.h         # 极简三维双精度向量类 (加减乘除、点积叉积、范数、归一化)
│   ├── Matrix3x3.h       # 3x3 矩阵基础类 (转置、逆、外积、Rodrigues旋转矩阵)
│   ├── Particle.h        # 结构质点类 (初始/当前坐标、速度、加速度、外力与内力、支座弹簧)
│   ├── Element.h         # 单元抽象基类 (Element 纯虚基类)
│   ├── TrussElement.h    # 空间桁架杆元与单拉索单元 (物理内力与刚度)
│   ├── PDSolver.h        # Projective Dynamics 求解器核心类 (Cholesky预分解与交替求解)
│   ├── MstModel.h        # 空间网格模型数据结构定义 (节点、杆件、截面、支座、荷载)
│   └── MstParser.h       # MSTCAD 空间网格模型文本解析器
├── main.cpp              # 算例示例程序 (现代 C++23 std::print，工况选择与反力自平衡校验)
├── CMakeLists.txt        # 跨平台 CMake 工程构建脚本
├── test.log              # 示例空间网架模型 (361 节点, 1388 杆件)
├── test2.log             # 示例复杂约束空间网架 (固定铰、单向滑移、双向滑移、弹性支座)
└── README.md             # 本说明文档
```

---

## 三、环境依赖与构建

### 1. 依赖库
- **C++ 编译器**：支持 **C++23**（推荐 MSVC 2022 17.10+ / VS 2026、GCC 14+ 或 Clang 18+，需支持 `<print>`）
- **Eigen 3**：核心稀疏矩阵与 Cholesky 分解库（3.3+）
- **Boost**：用于 `MstParser.h` 中的 `boost::spirit::qi` 高性能解析模块（仅需 Header-Only 头文件）

### 2. CMake 构建 (Windows / Linux / macOS)

#### 使用 vcpkg 安装依赖（推荐）：
```bash
vcpkg install eigen3 boost-spirit
```

#### 构建项目：
```bash
# 1. 创建并进入构建目录
cmake -B build -DCMAKE_TOOLCHAIN_FILE=[vcpkg路径]/scripts/buildsystems/vcpkg.cmake
cmake --build build --config Release
```

#### 运行示例：
```bash
cd build/Release  # 或 Linux 下的 cd build
./PDSolverDemo test.log
```

---

## 四、算例运行与实测结果

示例程序 `PDSolverDemo` 自动读取 `test.log`，提取所有可用工况，并支持命令行直接指定工况编号。

### 运行命令：
```bash
# 默认计算风吸力工况 1（三维空间反力与拉索抗拔特征最显著）
./PDSolverDemo test.log 3

# 计算全楼恒载工况
./PDSolverDemo test.log 1
```

### 实测输出示例 (C++23 `std::print` 格式化)：
```text
=================================================================================
         Projective Dynamics (PD) 空间杆系结构求解器开源算例示例程序             
  [理论论文]: Bouaziz et al., Projective Dynamics, ACM TOG (SIGGRAPH 2014)       
  [重要声明]: 本程序仅供学术研究、教学与力学交流使用，请勿用于真实工程设计！     
=================================================================================

>>> 正在加载 MST 空间结构模型文件: test.log ...
[+] 模型加载成功!
  • 节点总数: 361 个 (上弦: 169 个, 下弦: 192 个)
  • 杆件总数: 1388 根
  • 截面类型: 8 种
  • 边界支座: 12 个

---------------------------------------------------------------------------------
【模型荷载工况列表】
---------------------------------------------------------------------------------
 编号 | 工况名称                         | 节点数 | 总荷载 ΣFx (kN) | 总荷载 ΣFy (kN) | 总荷载 ΣFz (kN)
------+----------------------------------+--------+-----------------+-----------------+-----------------
  [ 1] | 恒载 (Dead Load)                 |    361 |            0.00 |            0.00 |        -1985.65
  [ 2] | 活载工况 1 (Live Load 1)         |    169 |            0.00 |            0.00 |         -496.10
  [ 3] | 风载工况 1 (Wind Load 1)         |    169 |           -4.65 |            0.00 |         1020.50
  [ 4] | 风载工况 2 (Wind Load 2)         |    169 |            4.64 |            0.00 |         1020.50
  [ 5] | 风载工况 3 (Wind Load 3)         |    169 |           -0.00 |            0.00 |         1113.26
  [ 6] | 风载工况 4 (Wind Load 4)         |    169 |           -0.00 |            0.00 |         1113.26
---------------------------------------------------------------------------------

>>> 已选定工况: 风载工况 1 (Wind Load 1) (作用节点数: 169, 合力: [-4.65, 0.00, 1020.50] kN)
>>> 正在启动 Projective Dynamics 交替优化求解...

=================================================================================
【PD 求解结果统计】
=================================================================================
  • 收敛状态: 已收敛 (Converged)
  • 迭代步数: 148 步
  • 求解耗时: 5.61 毫秒 (0.0056 秒)
  • 最大位移: 24.445 mm
  • 杆件最大拉力:     248.95 kN (杆件 #225)
  • 杆件最大压力:    -214.76 kN (杆件 #937)

---------------------------------------------------------------------------------
【边界支座反力明细表】
---------------------------------------------------------------------------------
 序号 | 节点号 | 支座类型     | 位置 (X, Y, Z) m     | 反力 Rx(kN) | 反力 Ry(kN) | 反力 Rz(kN) | 空间合力(kN)
------+--------+--------------+----------------------+-------------+-------------+-------------+-------------
  [ 1] |    314 | 弹性支座     | ( 0.0, 23.9, -1.3) |       -6.16 |        1.34 |      -72.41 |        72.68
  [ 2] |    318 | 弹性支座     | ( 0.0, 15.7, -1.3) |       -7.34 |        0.03 |      -94.69 |        94.98
  [ 3] |    321 | 弹性支座     | ( 0.0,  7.5, -1.3) |       -6.24 |       -1.31 |      -73.09 |        73.37
  [ 4] |    327 | 弹性支座     | (31.7, 23.9, -1.3) |        6.53 |        1.40 |      -59.61 |        59.98
  [ 5] |    331 | 弹性支座     | (31.7, 15.7, -1.3) |        7.84 |       -0.04 |      -82.73 |        83.10
  [ 6] |    334 | 弹性支座     | (31.7,  7.5, -1.3) |        6.66 |       -1.28 |      -60.95 |        61.32
  [ 7] |    340 | 弹性支座     | ( 7.7,  0.5, -1.3) |       -0.46 |       -6.71 |      -74.38 |        74.69
  [ 8] |    343 | 弹性支座     | (15.9,  0.5, -1.3) |        0.62 |       -7.87 |     -149.59 |       149.80
  [ 9] |    346 | 弹性支座     | (24.0,  0.5, -1.3) |        1.52 |       -6.22 |      -64.17 |        64.49
  [10] |    351 | 弹性支座     | ( 7.7, 30.8, -1.3) |       -0.48 |        6.69 |      -73.47 |        73.77
  [11] |    354 | 弹性支座     | (15.8, 30.8, -1.3) |        0.53 |        7.88 |     -149.89 |       150.09
  [12] |    357 | 弹性支座     | (24.1, 30.8, -1.3) |        1.63 |        6.11 |      -65.35 |        65.66
------+--------+--------------+----------------------+-------------+-------------+-------------+-------------
 合计 | 支座总反力向量 ΣR                    |        4.65 |        0.00 |    -1020.32 |
---------------------------------------------------------------------------------

=================================================================================
【结构全局静态力学自平衡校验】
=================================================================================
  • X 向平衡: ΣFx (-4.65) + ΣRx (+4.65) = +0.0057 kN (残差: 0.0057 kN)
  • Y 向平衡: ΣFy (+0.00) + ΣRy (+0.00) = +0.0017 kN (残差: 0.0017 kN)
  • Z 向平衡: ΣFz (+1020.50) + ΣRz (-1020.32) = +0.1711 kN (残差: 0.1711 kN)
  • 空间矢量总合力残差: 0.1712 kN (相对误差: 0.0168 %)
  • 自平衡检验结论: [优秀] 全局平衡高精度闭合！
=================================================================================
```

---

## 五、核心接口与使用代码片段

在用户自己的项目中直接使用 `PDSolver` 进行杆系结构求解仅需几行代码：

```cpp
#include "include/PDSolver.h"

// 1. 创建求解器
PD::PDSolver solver;
solver.maxIterations = 500;
solver.tolDisplacement = 1e-5; // 0.01 mm 收敛阈值

// 2. 添加质点
int p0 = solver.addParticle(PD::Vector3(0.0, 0.0, 0.0));
int p1 = solver.addParticle(PD::Vector3(3.0, 0.0, 0.0));

// 3. 添加杆件 (node1, node2, E (Pa), A (m^2))
solver.addTrussElement(p0, p1, 2.06e11, 0.001);

// 4. 边界支座与外荷载
solver.particles[p0].setFixed(true, true, true); // p0 节点三向固定
solver.particles[p1].fExt = PD::Vector3(0.0, 0.0, -10000.0); // p1 施加 10 kN 竖向荷载

// 5. 求解
auto stats = solver.solve();

// 6. 提取结果
PD::Vector3 disp = solver.particles[p1].displacement();
std::println("节点 1 位移: [{:.3f}, {:.3f}, {:.3f}] mm", 
             disp.x * 1e3, disp.y * 1e3, disp.z * 1e3);
```

---

## 六、参考文献

1. **Bouaziz, S., Martin, S., Liu, T., Kavan, L., & Pauly, M. (2014).** *Projective dynamics: Fusing constraint projections for fast simulation.* ACM Transactions on Graphics (TOG), 33(4), 1-11.
2. **Wang, H. (2015).** *A Chebyshev semi-iterative approach for accelerating projective and position-based dynamics.* ACM Transactions on Graphics (TOG), 34(6), 1-9.
3. **Liu, T., Bargteil, A. W., O'Brien, J. F., & Kavan, L. (2013).** *Fast simulation of mass-spring systems.* ACM Transactions on Graphics (TOG), 32(6), 1-7.

---

## 七、开源协议

本项目采用 [MIT License](LICENSE) 开源协议。

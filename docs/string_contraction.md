# 从 unique-spin-string 直接收缩到 VB structure

> 算法整理与共同审阅稿 · v1 · 2026-10-07
>
> 本文讨论固定轨道点上的精确收缩，以及这一收缩的反向微分。先把表示和求和顺序讲清楚，不把 RI、正交 direct-CI、轨道优化器等其他贡献混进来。

## 1. 我们究竟改变了什么？

VB 结构由成键模式和自旋耦合定义。一个结构通常包含许多完整行列式展开项；不同结构的展开项，又大量复用相同的单自旋占据串。

传统计算可以先计算这些单自旋核，再逐个组合出完整行列式对的矩阵元，最后根据展开系数累加到结构矩阵。我们的核心改变是：**把结构展开系数组织为 spin-string 系数块，直接在这些块之间收缩，不再在主 H/S 收缩与 action 中逐个形成并散射完整行列式对的标量。**

这改变的是精确计算的组织方式，不改变所选 VB 结构、结构系数的变分空间或给定积分核。

```mermaid
flowchart LR
    K[共享的单自旋算子核]
    K --> O[旧路径：组合完整行列式对标量]
    O --> P[按归属与展开符号散射]
    P --> S[VB structure 矩阵或向量]
    K --> N[新路径：与结构系数块直接收缩]
    N --> S
```

这里有一个重要的历史边界：XMVB4.0 已经复用 unique-spin-string 核，并非对每个完整行列式对都从头计算单自旋矩阵元；它也没有在本文核对的普通 RDM–VBSCF 路径中保存完整的行列式 H/S 数值矩阵。我们优化的是核之后的组合与结构散射层。

## 2. 三个层级，不能混为一谈

| 层级 | 表示什么 | 在本文算法中的作用 |
| --- | --- | --- |
| VB structure | 指定的成键／自旋耦合波函数构件 | 真正的结构系数变分空间 |
| 完整 determinant | 一条 α 串与一条 β 串组成的展开项 | 描述结构展开的组合关系，不必成为数值矩阵元中间层 |
| unique-spin-string | 去重后的单自旋占据串 | 共享单自旋算子核的索引空间 |

在固定、非正交轨道下，记 α、β 串分别为 $a$、$b$，结构编号为 $I$。采用实数系数约定：

$$
|\Phi_I\rangle
=\sum_{a,b}(C_I)_{ab}|a_\alpha b_\beta\rangle,
\qquad
|\Psi\rangle=\sum_I v_I|\Phi_I\rangle.
$$

$C_I$ 是结构展开的精确系数矩阵。它不是新的变分参数，也不是对结构波函数的近似。

因此，结构向量 $v$ 对应的 string-product 系数为

$$
X(v)=\sum_I v_I C_I.
$$

内部出现 $X$，不意味着优化器在整个 string-product 空间任意改变波函数。每次 action 都必须投影回原来选定的结构空间。

## 3. 一个贯穿全文的例子：两对 singlet

### 3.1 结构展开与费米符号

考虑四个轨道上的两对 singlet。定义未归一化配对产生算符

$$
P_{ij}^{\dagger}
=a_{i\alpha}^{\dagger}a_{j\beta}^{\dagger}
-a_{i\beta}^{\dagger}a_{j\alpha}^{\dagger}.
$$

选择结构

$$
|\Phi_A\rangle=-P_{12}^{\dagger}P_{34}^{\dagger}|0\rangle.
$$

约定完整行列式的产生算符按“α 在前、β 在后；各自轨道编号升序”排列。上式的整体负号只用于选择统一的结构相位；不能忽略重排产生算符产生的费米符号。

| 展开项 | α 占据串 | β 占据串 | 系数 |
| --- | --- | --- | ---: |
| $D_1$ | $(1,3)$ | $(2,4)$ | $+1$ |
| $D_2$ | $(1,4)$ | $(2,3)$ | $+1$ |
| $D_3$ | $(2,3)$ | $(1,4)$ | $+1$ |
| $D_4$ | $(2,4)$ | $(1,3)$ | $+1$ |

这是 **4 个完整展开项、4 条 unique α strings、4 条 unique β strings**，不是两条 α、两条 β。

### 3.2 展开系数变成一个稀疏连接图

按如下局部顺序排列：

$$
\alpha_A=(13,14,23,24),\qquad
\beta_A=(24,23,14,13).
$$

那么

$$
C_A=\begin{pmatrix}
1&0&0&0\\
0&1&0&0\\
0&0&1&0\\
0&0&0&1
\end{pmatrix}.
$$

这里的对角形式来自人为选择的匹配列序，不代表一般结构的系数矩阵必然对角。实际存储可以只保存非零连接：

```text
α13 ──(+1)── β24
α14 ──(+1)── β23
α23 ──(+1)── β14
α24 ──(+1)── β13
```

第二个结构可选

$$
|\Phi_B\rangle=-P_{13}^{\dagger}P_{24}^{\dagger}|0\rangle.
$$

其展开是：

| α 占据串 | β 占据串 | 系数 |
| --- | --- | ---: |
| $(1,2)$ | $(3,4)$ | $+1$ |
| $(1,4)$ | $(2,3)$ | $-1$ |
| $(2,3)$ | $(1,4)$ | $-1$ |
| $(3,4)$ | $(1,2)$ | $+1$ |

不同结构的局部 string 列表可能不同。计算交叉矩阵元时，必须按照左右 support 提取对应的核，或者把两者嵌入统一的全局索引。不能因为两个局部矩阵尺寸相同，就把相同行号当成相同的 string。

## 4. 重叠：从四重求和到直接结构收缩

### 4.1 原始表达式

定义未归一化单自旋重叠核

$$
S^\alpha_{aa'}=\langle a_\alpha|a'_\alpha\rangle,
\qquad
S^\beta_{bb'}=\langle b_\beta|b'_\beta\rangle.
$$

完整展开项之间的重叠为两者的乘积。因此

$$
S_{IJ}
=\sum_{a,a',b,b'}
(C_I)_{ab}S^\alpha_{aa'}(C_J)_{a'b'}S^\beta_{bb'}.
$$

旧的组织方式可以逐个构造
$S^\alpha_{aa'}S^\beta_{bb'}$，再乘结构系数、累加。两对 singlet 的结构自重叠有 16 个有序展开项对；实际实现可以进一步利用对称性。

### 4.2 精确重排求和

先对右侧 β 指标求和，再处理 α 指标：

$$
W=C_J(S^\beta)^T,\qquad
Y=S^\alpha W.
$$

用实数 Frobenius 内积

$$
\langle L,R\rangle_F=\operatorname{Tr}(L^T R)
=\sum_{a,b}L_{ab}R_{ab},
$$

就得到核心公式

$$
\boxed{S_{IJ}=\langle C_I,S^\alpha C_J(S^\beta)^T\rangle_F.}
$$

| 中间量 | 左右结构局部 support 下的形状 |
| --- | --- |
| $C_J$ | $U_{\alpha,J}\times U_{\beta,J}$ |
| $S^\alpha$ | $U_{\alpha,I}\times U_{\alpha,J}$ |
| $S^\beta$ | $U_{\beta,I}\times U_{\beta,J}$ |
| $W$ | $U_{\alpha,J}\times U_{\beta,I}$ |
| $Y$ | $U_{\alpha,I}\times U_{\beta,I}$ |

这一等式没有做截断。优化来自求和次序、稀疏连接与块计算的利用，而不是来自舍弃展开项。

## 5. Hamiltonian 仍使用同一套收缩结构

### 5.1 一电子与 same-spin 双电子通道

记 $H^\alpha$、$H^\beta$ 包含各自的有效一电子及同自旋双电子贡献，则

$$
H^{\mathrm{same}}_{IJ}
=\left\langle C_I,
H^\alpha C_J(S^\beta)^T
+S^\alpha C_J(H^\beta)^T
\right\rangle_F.
$$

不活性空间的标量能量另贡献 $E_{\mathrm{core}}S_{IJ}$；不能将它遗漏，也不能在有效核和标量项之间重复计算。

### 5.2 Opposite-spin 通道

定义未除以 string overlap 的单自旋跃迁密度核

$$
D^{\sigma,ij}_{aa'}
=\langle a_\sigma|a_{i\sigma}^{\dagger}a_{j\sigma}|a'_\sigma\rangle.
$$

将轨道对 $(i,j)$ 简记为通道 $p$，并将积分及所选索引约定包含在 $V_{pq}$ 中，则

$$
H^{\mathrm{opp}}_{IJ}
=\sum_{p,q}V_{pq}
\left\langle C_I,
D^{\alpha,p}C_J(D^{\beta,q})^T
\right\rangle_F.
$$

因此，一条通用通道都具有

$$
F_{IJ}=\langle C_I,A C_J B^T\rangle_F
$$

的形式。Exact 与 RI 会改变通道如何表示、组合和应用，但不会改变两端的结构投影关系。本文不据此声称所有通道都具有相同的计算成本。

## 6. Davidson：在结构空间求解，在 string 空间执行 action

真正要求解的广义本征问题是

$$
H_{\mathrm{str}}v=E S_{\mathrm{str}}v,
\qquad v\in\mathbb R^{N_{\mathrm{str}}}.
$$

不必先组装完整结构 H/S，就可以对结构向量执行：

```mermaid
flowchart LR
    V[结构向量 v] --> X[展开为 string 系数 X]
    X --> A[单自旋通道与 X 收缩]
    A --> Y[string-space image]
    Y --> P[与各结构系数块投影]
    P --> R[结构向量 Hv 与 Sv]
```

对 overlap 而言：

$$
X=\sum_J v_J C_J,\qquad
Y_S=S^\alpha X(S^\beta)^T,\qquad
(Sv)_I=\langle C_I,Y_S\rangle_F.
$$

Hamiltonian 按上述物理通道累加，得到 $Y_H$，再计算

$$
(Hv)_I=\langle C_I,Y_H\rangle_F.
$$

这种 expand–action–project 精确表示所选结构子空间上的算子。内部 image 不必仍是某个单独结构，也不需要所选结构在每个内部操作下闭合；关键是输入映射和最终投影始终使用原来的 $C_I$。

概念伪码如下。它说明线性映射，不要求实现分配完整的稠密中间量：

```text
apply_structure_block(V):
    X = expand_selected_structures(V)
    Y_H, Y_S = zero_images()

    for each physical channel:
        accumulate_string_action(X, Y_H, Y_S)

    return project_to_selected_structures(Y_H, Y_S)
```

| 核对的模块 | 求解空间 | 方法 |
| --- | --- | --- |
| XMVB4.0 普通 RDM–VBSCF | VB structure | dense 广义对角化 |
| XMVB4.0 FCI 模块 | 完整 determinant | Davidson |
| 本项目 VBSCF Davidson | 所选 VB structure | matrix-free Davidson |

FCI 使用行列式空间是其问题定义决定的，并不是错误。不能把官方 FCI 模块与普通 selected-VB 求解混为一谈。

## 7. Backward：直接生成单自旋核的反向权重

考虑实数通用收缩

$$
F=\langle L,A R B^T\rangle_F.
$$

固定 $L,R$，其微分为

$$
dF=\langle L,(dA)R B^T+A R(dB)^T\rangle_F.
$$

通过迹的循环性质得到

$$
\boxed{\overline A=\frac{\partial F}{\partial A}=LBR^T,\qquad
\overline B=\frac{\partial F}{\partial B}=L^TAR.}
$$

也就是说，结构层提供的反向权重可以通过系数块直接收缩到 α、β 核。不必先生成每个完整行列式对的导数，再依展开关系散射。

对于归一化、驻定的单态结构系数 $v$，固定点能量导数包含

$$
dE=v^T(dH-E\,dS)v.
$$

态平均则按态权重组合相应贡献。这里仍需正确处理归一化与结构驻定条件；两条通用 adjoint 公式不是完整轨道梯度的全部内容。

HVP 可以继续对收缩求方向微分。例如

$$
\dot{\overline A}
=\dot LBR^T+L\dot B R^T+LB\dot R^T.
$$

核的方向响应、结构系数响应以及轨道坐标的微分，需要在完整 HVP 中一并处理。因此，我们的主张是**收缩及其微分可以在 string 层闭合**，不是结构响应消失了。

## 8. Tile：计算组织和内存组织必须一起考虑

上述矩阵公式说明代数，不强制保存所有 string-pair 核。可以按核块生成和消费：

```text
for each spin-pair tile:
    generate the kernels required by this tile
    contract with the relevant coefficient support / input block
    accumulate structure images or kernel adjoints
    release or reuse tile workspace
```

Tile 的关键是让核在有效生命周期内完成所有必要消费，避免为了省内存而在多个消费者中重复生成。具体通道可能需要不同的 sweep 和 workspace；这个伪码不是声称所有消费者已经融合为一次遍历。

| 对象 | 可行的组织方式 | 仍需警惕 |
| --- | --- | --- |
| 单自旋 pair 核 | 按 tile 生成／消费 | 全局 accepted 或 directional cache |
| 单个结构的展开系数 | 稀疏连接与局部 support | 初始化的展开项数量 |
| Davidson 输入与 image | block action、分块投影 | 完整 spin-product carrier 的尺寸 |
| 结构 H/S | matrix-free action | dense 对照路径仍需要结构矩阵 |

Tile 能限制局部 pair 载荷，但不自动证明整个程序内存与 $U_\alpha,U_\beta$ 无关。

## 9. 计算标度：哪里降低了，哪里没有？

### 9.1 可分离 action 的一个清晰比较

对完整稠密 spin-product 输入，单通道逐 pair 应用可有

$$
O(U_\alpha^2U_\beta^2)
$$

成本。将同一通道写成 $A X B^T$ 后，稠密矩阵计算成本为

$$
O(U_\alpha^2U_\beta+U_\alpha U_\beta^2).
$$

若 $U_\alpha\sim U_\beta\sim U$，这一特定 action 从四次变为三次复杂度。Davidson block width 为 $b$ 时，上述算术量相应乘以 $b$。

这个比较不包括核生成、通道数量、积分应用及完整梯度/HVP 成本。也不能把完整稠密输入的比较直接套到一个稀疏单结构自矩阵元上。

### 9.2 为什么不能说“指数墙已经消失”？

对由 $n$ 对不相交 singlet 构成的单结构，一般有 $2^n$ 个展开项；这种例子本身也有 $2^n$ 条 α strings 和 $2^n$ 条 β strings。

新算法避免的是把这些连接在数值热路径中逐对组合、查询、散射，并利用跨结构的共享核及矩阵收缩。**它并未证明 unique-string 数量或初始化展开变成多项式。**

因此更准确的方法主张是：

> 在保持所选 VB 结构空间精确性的前提下，将主要数值收缩从完整行列式对的组合与散射，重组为共享 spin-string 核上的直接结构收缩；由此支持高效块计算、反向收缩及结构空间 matrix-free 求解。

是否带来数量级总加速，需要隔离 kernel 和端到端计算两类证据，不能只凭某个理想 action 的 FLOP 估算。

## 10. 已核对的教学例子

对第 3 节的结构 A、B，取正定、非正交的单电子重叠矩阵

$$
O=\begin{pmatrix}
1&0.18&0.07&0.03\\
0.18&1&0.11&0.06\\
0.07&0.11&1&0.16\\
0.03&0.06&0.16&1
\end{pmatrix}.
$$

由相应占据子矩阵的行列式构造单自旋重叠核，分别采用完整展开求和与矩阵收缩。

| 检查 | 结果 |
| --- | ---: |
| 16 个项对直接求和的 $S_{AB}$ | $-2.0759798600$ |
| string 矩阵收缩的 $S_{AB}$ | $-2.0759798600$ |
| 两种算法差值 | $0$（本次双精度计算） |
| α 核 adjoint 的方向有限差分绝对误差 | $1.24\times10^{-9}$ |
| β 核 adjoint 的方向有限差分绝对误差 | $1.54\times10^{-9}$ |

这里计算的是未归一化结构的交叉重叠，负值不表示错误。它是教学代数检查，不是实际分子性能测试，也不替代完整梯度/HVP 验证。

检查程序目前位于 [生成与检查脚本](../tools/render_string_contraction.py)，原始结果见 [检查记录](string_contraction_checks.json)。该脚本同时生成上一版分页 PDF；本文是后续审阅的主要文档。

## 11. 代码证据与下一轮要共同确认的问题

以下行号固定到审阅时的 C++ 快照 `d25d96c4b9d281328e68212f6ffa83052e23e93c`；后续代码改变后应重新核对。官方来源为附件仓库的纯净 `xmvb4.0` 分支 `a3613ed32ebffcaf2d98c1c9124555708bea995b`，不是 DeepVBH。

| 命题 | 文件／位置 |
| --- | --- |
| 官方单自旋核、pair 组合、结构散射 | `src/vbscf/gradient_rdm.F90`：`Hamhd`、`Hov1`、`Hovstr1`，约 2111、2876、3380 行 |
| 官方普通 VBSCF dense 广义对角化 | `src/math/eigencalc.c:85`，`LAPACKE_dsygv` |
| 官方 FCI Davidson 使用 `det_num` | `src/mc_scf/fci_cal.c:467` |
| 当前直接系数块收缩 | `src/vbscf/structures/assembly/block_kernels.hpp:31` |
| 当前结构 expand／project | `src/vbscf/structures/assembly/action.cpp:858`／`:778` |
| 当前初始化仍展开配对结构 | `src/vbscf/structures/expansion/expander.cpp:66` |

主路径的结论不能扩展成“代码任何地方都没有完整项对计算”：例如某些对角／预条件器设置仍可能按需访问项对。

下一轮建议先讨论：

1. 第 3–4 节是否已经让读者理解三个层级、费米符号和求和重排？
2. 我们最希望强调的是跨结构共享、可分离 action，还是 backward 的直接收缩？哪些应分别给出证据？
3. 是否补一个两个结构共享同一条 string 的图，将“跨结构复用”画得更直观？
4. 如何设计同输入、同精度的隔离测速，使结果真正归因于收缩算法，而不是积分技术或优化器？

这些问题留给共同审阅，不先用未经验证的速度数字替代推导。

# 2C-VBSCF implementation plan

理论依据：[Two-Component VBSCF Theory](docs/jctc_2c_vbscf.tex)。

状态：**计划阶段，尚未实现两分量计算。** 不将现有非正交 VBSCF、计算正交化或 scalar-relativistic 修正称为 2C-VBSCF。

命名：通用方法为 **2C-VBSCF**；采用明确的 X2C 有效哈密顿量时为 **X2C-VBSCF**。

## 目标和基本原则

- 保留价键结构的局域占据、配对和自旋耦合母结构；波函数仍然是所选 VB 结构的线性组合。
- 采用共享的非正交复数 spinor，支持 HAO 严格 AO 支持域与 OEO full-AO 支持域。
- 通过 occupation strings、因子化 exterior image 和精确结构投影计算 H/S action，不建立完整的行列式对矩阵元计算层。
- 同时支持 exact Coulomb 与 RI；RI 不是两分量方法定义的一部分。
- 尽量复用已经验证的代数和优化流程，但不将整个代码库机械替换成复数。
- 命名短而明确，Eigen 表示物理矩阵／向量，Doxygen 解释维数、打包、共轭与导数约定。
- 不添加自动切回 spin-free、cofactor 或 dense 的 fallback。未支持的模型、超出内存预算的 carrier、奇异 frame 和不完整输入应明确报错。
- 每一阶段独立验收、提交 commit；不保留临时生产分支、无效参数或重复版本。
- 构建和数值计算在 **瀚海25的 Slurm 计算节点**执行。这里只做源码审查和文档编辑；本机及登录节点不运行计算任务。

## 当前代码中的实际边界

| 当前模块 | 可以复用 | 必须推广或重新验证 |
| --- | --- | --- |
| `src/vbscf/structures/expansion/types.hpp` | 结构标签、占据拓扑、稀疏结构系数映射的思路 | 现有数据固定单一 multiplicity 与 alpha/beta 电子数；不能直接表达全部 spinor 扇区和耦合母结构 |
| `src/vbscf/structures/orthogonal_ci/exterior_transform.*` | 因子化 exterior 操作、连接符号、方向导数思路 | 实数、单个 fixed-spin carrier、普通转置；通用 spinor 变换还会跨计算分组 |
| `src/vbscf/structures/orthogonal_ci/sigma.*` | string link、block action、积分 adjoint 的组织 | spin-free alpha/beta 分解、实数 pair kernel、实数积分对称性 |
| `src/vbscf/structures/orthogonal_ci/integrals.*` | Gram 分解与 moving-frame 求导链 | 复数 Hermitian 分解、ordered-pair 积分、两分量有效算子 |
| `src/core/eigensolver.*` | Davidson 增量子空间、H/S image 复用、残差认证 | 当前向量、action 和内积为实数；需要复数 Hermitian 版本 |
| `src/core/eigen_response.*` | bordered response、等权 subspace response、同点 recycle | Hermitian 投影、复数规范、Kramers/简并簇；不能直接使用实数转置 |
| `src/vbscf/orbitals/charts/*` | 支持域与局部坐标的设计原则 | 不能原样继承单轨道归一化、相位、inactive/active 冗余结论 |
| `src/vbscf/optimization/driver/*` | 接受点／试探点流程、停止条件、line search/NEO 组织 | 必须通过已验证的两分量 objective、gradient、HVP、metric 接口接入 |

现有 scalar VBSCF 是一个独立、仍需维护的物理模型，不是两分量计算失败时的兼容路径。本工作不同时进行无关目录重构或删除 scalar 方法。

## 模型与数值契约：先确定，再写内核

### 波函数

- 每个空间局域标签对应两个 spinor 构件；构件标签不代表物理 alpha/beta。
- 结构字典明确记录电荷／配对结构与耦合母结构，不只记录一个总 multiplicity。
- 支持任意显式选择的结构字典；不自动补齐或裁剪用户的变分空间。
- 完整自旋通道的构造与进一步 selected-channel 截断分开表示。多重度、投影分量和耦合路径均由字典生成规则明确指定。
- 默认理论覆盖独立 spinor；Kramers-paired 是一个明确的约束版本。第一轮同时用微型 fixtures 检验二者，不用约束模式冒充一般模型。
- 普通 VBSCF 回归必须限制到同一 collinear 波函数及同一结构字典。仅关闭 SOC、仍允许一般 spinor，并不保证变分空间相同。

### 数值表示

- AO overlap 保留实数；spinor 系数、变换、结构振幅、H/S image 和响应使用 `Eigen::MatrixXcd`/`VectorXcd`。
- 本征值、Hermitian 对角、残差范数与优化器坐标使用实数。
- 复数内积明确使用共轭。定义实标量导数与 complex adjoint 的配对约定，统一检查实部和虚部方向。
- 计算分组使用 A/B 标签，而非误称 alpha/beta；固定电子数内允许全部所需分组扇区。
- 全局轨道排序、bitstring 位序、string 编号、矩阵 leading dimension 及 block 打包只有一种定义。
- 局部维数命名使用 `n_bf`（空间 AO 基函数）、`n_orb`（局域轨道标签）、`n_spinor`、`n_active_orbitals`、`n_active_electrons`；不以含糊的 `n_ao` 同时表示 AO 与 active orbitals。
- 不把一般复数双电子积分压入实数 symmetric-pair/八重对称格式。
- spinor frame 与结构空间分别检查线性独立性。零 occupied-pair overlap 不等同于完整 frame 奇异。
- 原则上不以正则化或 pseudoinverse 掩盖未知的结构／frame 秩变化。

## 分阶段里程碑

### P0 — 独立微型参考与契约

- [ ] 定义上述数据契约、字典及索引约定。
- [ ] 建立独立稠密 fermion reference，仅用于小空间测试，不进入生产路径。
- [ ] 建立 `(2o,2e)` 六维完整空间及任意三个结构的 fixtures。
- [ ] fixtures 包含随机复数 h、符合 Hermitian/交换对称性的 Coulomb 张量、可逆非正交 spinor frame、稀疏结构映射，以及可控的零 occupied-subset overlap。
- [ ] 增加一个奇电子微型模型与一个等权简并子空间模型，检验时间反演与子空间规范。
- [ ] reference 用独立的显式创建／湮灭代数，不调用被测的 string link 或 exterior 内核。

验收：fixtures 自身通过 Hermitian 性、粒子数、独立结构秩和已知解析极限检查；提供模型、种子、维数、Hamiltonian conventions 的记录。

### P1 — 通用 occupation-string 算子

- [ ] 构造 A/B 多扇区 strings 与跨组创建／湮灭连接。
- [ ] 实现复数 one-body block action，统一处理普通和 spin-mixing 矩阵元。
- [ ] 实现 exact two-body action 与正确的 normal-ordering 减项。
- [ ] 实现 RI 两次 one-body action 与 one-body counterterm，共享连接及 workspace。
- [ ] 仅在确有共同代数时提取共享连接逻辑，不复制完整旧 sigma 再大量加条件分支。

验收：所有列及随机 block 与 P0 reference 一致；覆盖跨组一步／两步转移、重复指标、对角项、Hermitian adjoint、不同 block width。RI 使用同一批 factors 重建的张量作为 reference，不能拿拟合误差判断内核错误。

### P2 — 非正交 exterior image 与 VB 精确投影

- [ ] 实现复数 Gram/Cholesky computational frame。
- [ ] 推广因子化 exterior transform，包含跨 A/B 分组的变换；不能只将旧左右两个独立 transform 改成复数。
- [ ] 实现结构 `expand`、image transform、Hamiltonian action、conjugate-adjoint contraction。
- [ ] 实现 H/S block action；任意所选结构的整个 image 均需被计算 carrier 容纳。
- [ ] 支持固定 spinor frame 下的 Kramers-paired fixtures，并检查结构空间的时间反演闭包。

验收：完整及任意三个结构的 H/S action 与独立非正交 reference 一致；投影前后 Rayleigh 商一致；不得依靠 selected-space 闭包假设。一般变换应覆盖原始结构支持以外的占据分量。

### P3 — 复数 Hermitian 结构求解器

- [ ] 为 core generalized dense/Davidson 提供复数数值契约。
- [ ] 保留一个共享的已验证求解流程；标量类型推广限定在必要的线性代数层，不建设泛化框架。
- [ ] 保留 H/S image 的增量更新、metric orthogonalization、warm start 和 residual certification。
- [ ] 结构 diagonal 用实数 Hermitian 对角，不取错误的逐元素绝对值替代。
- [ ] 增加多根及等权简并簇的检查；比较子空间 projector，不比较任意相位的单列向量。

验收：dense/Davidson 的低根能量、metric orthonormality、CI stationarity 与 P0/P2 reference 一致；测试 H 加常数倍 S 后残差判断的稳健性。生产 Davidson 不构造完整 H/S，不因最终输出而强制全谱求解。

### P4 — 分子静态 X2C Hamiltonian 与输入定义

- [ ] 定义 AO-level `h_2c` 提供接口与基组顺序、归一化、分量排序、能量零点的 provenance。
- [ ] 先从可信参考生成 Hamiltonian fixtures；导入仅是 reference 路径，不是完成生产 X2C 的标志。
- [ ] 实现并验证所选原生 X2C AO 算子提供方式。是否使用新的依赖需先核查许可、构建及集群可用性。
- [ ] 将 X2C1e 与指定的静态 AMF 模型区分，不混加另一份 SOC。第一版不隐式引入 density-dependent molecular SOMF。
- [ ] 复用 AO bare-Coulomb exact/RI providers，新增正确的复数 spinor 变换；不重新复制两套 AO ERI 存储。
- [ ] 明确 inactive/frozen-core spinor、core energy、有效 active h，以及 core orbital 导数。不能直接套用标量 doubly-occupied 公式。
- [ ] 定义实际输入语法与复数初猜格式，并进行 round-trip 检查；显式请求但缺失的初猜不得被静默替换。
- [ ] 明确 `(mo,Ne)` 与 `2m` 个 active spinor 的数量关系，输出同时报告两者。

验收：相同固定 AO 基、Hamiltonian approximation、core 和 spinor 下，对照独立 2C-CASCI/reference 的能量；把 Hamiltonian 差异与结构截断差异分开。普通 scalar PySCF CASSCF 不是一般 spinor 模型的 reference。

### P5 — 解析梯度与受约束坐标

- [ ] 实现 Gram、Cholesky、exterior image、h/exact/RI 变换的解析 directional derivative 与 adjoint。
- [ ] 包含左右 image 及 moving computational-basis operator 的全部导数，不重复求同一 accepted-point 中间量。
- [ ] 实现由结构能量及 overlap 导数得到的解析轨道梯度；CI stationarity 必须达到导数验收所需精度。
- [ ] 建立 real/imag 坐标和 support-preserving 更新；明确 independent/Kramers 版本的 pullback。
- [ ] 根据结构闭包认证 gauge，不自动删除 active--active、自旋构件相位或局域 frame 方向。
- [ ] 仅在证明可补偿后使用 normalization/gauge fixing。selected dictionary 的单 spinor scaling 不作为无条件冗余。
- [ ] 等权态平均使用同一 roots、weights 和 Hamiltonian，对接受点、试探点及有限差分均保持一致。

验收：HAO/OEO、independent/Kramers、selected/full 的实部与虚部能量有限差分；多步长展示截断区与舍入区。检查 adjoint pairing，而不只比较最终收敛能量。

### P6 — 完整 HVP、结构响应及 coupled 接口

- [ ] 从已经验证的梯度链产生解析 HVP，包括参数化/retraction 的二阶项。
- [ ] 推广复数 generalized eigen response 与等权 subspace response。
- [ ] 区分 fixed-coefficient orbital action、coupling、structure response 和完整 relaxed HVP，名称与理论定义一致。
- [ ] 对 isolated roots 使用投影／bordered residual；对简并/Kramers 簇使用整个子空间，不能穿过零 gap 使用单态公式。
- [ ] 同点 response space 与已计算 H/S images 可复用；更新轨道点后不能无条件保留旧算子 image。
- [ ] 初始正确性版本求到独立 residual 认证精度。响应近似、延迟 coupling 与误差分配在完整算子验证后再评估。

验收：HVP 对梯度有限差分，实坐标 Hessian 的对称性，coupled/Schur 两种表达的一致性，独立小型 dense Hessian 对照。全部包含结构归一化与 overlap 导数项。生产不要求完整结构本征谱。

### P7 — 接入现有优化器、输出与分子回归

- [ ] 先通过已验证的 objective/gradient/metric 接入现有一阶驱动，再接 TNHVP/NEO 的完整二阶接口；不另写一整套 SCF 生命周期。
- [ ] one-electron block preconditioner 如要用于 spinor 坐标，必须先验证其 realified block 与 constrained pullback，而非复用旧实数 block。
- [ ] retraction、horizontal lift、gradient pullback 和 HVP 属于同一参数化；跨 chart 的 secant/history 传输需明确处理。
- [ ] 不更改现有 `ISCF=5/7/8/9` 的算法含义。两分量 Hamiltonian、轨道约束与优化器选择是独立输入字段。
- [ ] 输出明确标识 2C/X2C 模型、core、支持域、parent channels、Kramers restriction、roots/weights、CI residual 和物理梯度。
- [ ] 可输出各态的自旋期望与混合信息，但 SOC 下不把整数 multiplicity 或 spin purity 当作最终态必需条件。
- [ ] 保留复数 spinor restart；验证重新读入的第零步能量和结构空间不变。

验收顺序：

1. 微型完整／selected 模型端到端收敛与 dense reference 一致。
2. 零 SOC 且限制到相同 collinear ansatz，恢复现有 scalar VBSCF 的能量和导数。
3. 小型 full-space OEO 对照使用同一模型的 2C-CASSCF/reference；不保证不同优化器必然落入同一 basin，需同时检查 stationarity 和局部 curvature。
4. 增加一个奇电子模型的 Kramers/等权态平均检查。
5. 再进入过渡金属 HAO/OEO 和多态测试，不以 F2 单例作为充分验证。

## 建议的文件组织

按现有功能分区增加少量叶目录，不建立第二套 `2c-vbscf` 总目录或第二个 application driver。

```text
src/vbscf/
  orbitals/spinor/
    frame.hpp/.cpp         # Gram、计算 frame、Kramers 约束
    coordinates.hpp/.cpp   # 支持域、实坐标、已认证 gauge、pullback
  structures/spinor/
    basis.hpp/.cpp         # 母结构字典与结构映射
    strings.hpp/.cpp       # A/B 扇区、索引与连接
    transform.hpp/.cpp     # 因子化 exterior image、伴随与导数
    action.hpp/.cpp        # exact/RI string action 与结构投影
  integrals/x2c/
    hamiltonian.hpp/.cpp   # 明确的 AO effective-operator 提供方式
  integrals/spinor/
    transform.hpp/.cpp     # 一／二电子与 RI factors 的复数变换
  derivatives/spinor/
    gradient.hpp/.cpp
    hessian.hpp/.cpp
```

以上是建议的最小职责划分，不要求提前创建空文件。每个文件随对应里程碑出现，新增 `.cpp` 只在现有 sources manifest 中归属一次。

core eigensolver/response 的必要推广保留在 `src/core`，不复制到 spinor 文件夹。真正可复用的底层数学例程可以抽取；不得为了避免重复而强迫 spin-free 与 general-spinor 使用不适合其连接结构的同一个内核。

测试进入现有 `tests/vbscf`，在 `targets.cmake` 注册。独立稠密 reference 只供测试，不成为生产 fallback。

## 内存与性能约束

- 不建立 carrier 的 `D × D` Hamiltonian/overlap，不建立 spin-string-pair 的全局 `U²` 物理缓存。
- 不显式建立 exterior compound matrix；保留因子化变换。
- dense 结构对角化保留为显式选择；分配前检查 `n_str²` 成本。Davidson 默认走结构 action。
- full carrier 至少有 `D × block_width` 的成本，复数元素为 16 bytes；所有扇区求和后再估算，不能只检查最大的一个扇区。
- 初始 planner 包含 carrier、H/S images、Davidson basis、连接表、transform links、RI factors 和每线程 workspace 的同时存活峰值。
- 连接表也可能超过振幅内存：避免为所有 spinor pairs、auxiliaries、roots 或线程复制同一拓扑。
- RI auxiliary 与 complex integral transformation 使用有界 workspace；不能同时常驻全部 factors、全四指标积分及其全部 directional 副本。
- 并行优先作用于合法的 block/tile 收缩；检查 scatter 竞争及归约 workspace，禁止每线程复制整个 carrier。
- 内存不足时拒绝请求并报告估算项；不能通过截断 image、丢弃自旋扇区或切换另一模型来“继续计算”。
- 性能测试在正确性通过后进行；报告步数、CI/actions、orbital/response actions、wall time、峰值 RSS，以及 carrier 实际维数。没有经过测量前不承诺加速比。

## 第一轮实际实施范围

**P0 → P1 → P2**：先完成复数 string 内核与任意 selected-VB 的精确固定轨道 action。

这轮不同时启动真实分子 X2C、轨道优化或大活性空间性能调参。结束时必须能证明：给定同一个非正交 spinor frame 和同一个 VB 字典，H/S action 与独立参考一致。

每轮记录完成的复选项、commit、Slurm job、测试输出位置与尚未通过的条件；没有验收的数据不得标记为已完成。

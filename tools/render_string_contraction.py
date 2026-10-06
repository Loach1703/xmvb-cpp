#!/usr/bin/env python3
"""Render the reviewable Chinese group-meeting deck (vector PDF + overview).

This is document generation, not a molecular benchmark. The small numerical
checks validate the teaching example and the displayed contraction identities.
Requires Python 3, matplotlib, numpy, Pillow, and a CJK font.
Rebuild: uv run --python /usr/bin/python3 --with matplotlib --with numpy
         --with pillow tools/render_string_contraction.py
"""

import argparse
import io
import itertools
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from PIL import Image, ImageOps

ROOT = Path(__file__).resolve().parents[1]
INK, MUTED = "#182C40", "#5F7285"
BLUE, TEAL, OLD = "#2563A6", "#087F82", "#AD6738"
PALE, LINE, WHITE = "#F0F5FA", "#D7E1EB", "#FFFFFF"
FONT = None
PAGES = []

# Render mathematical indices with the math font, not unsupported CJK glyphs.
MATH_LABELS = {
    "Σₐᵦ": r"$\sum_{ab}$", "Σᵢ": r"$\sum_I$",
    "(Cᵢ)ₐᵦ": r"$(C_I)_{ab}$", "Yₐᵦ": r"$Y_{ab}$",
    "(Hv)ᵢ": r"$(Hv)_I$", "Sₐᵦ": r"$S_{AB}$",
    "(Sβ)ᵀ": r"$(S^\beta)^T$", "LᵀR": r"$L^TR$",
    "αₐ, βᵦ": r"$\alpha_a,\beta_b$", "Sᵢⱼ": r"$S_{IJ}$",
    "Cₐ": r"$C_A$", "Cᵢ": r"$C_I$", "Cⱼ": r"$C_J$",
    "vᵢ": r"$v_I$", "αᵢ": r"$\alpha_I$", "βᵢ": r"$\beta_I$",
    "αⱼ": r"$\alpha_J$", "βⱼ": r"$\beta_J$",
    "D₁": r"$D_1$", "D₂": r"$D_2$", "D₃": r"$D_3$", "D₄": r"$D_4$",
}


def text(ax, x, y, value, size=20, color=INK, **kwargs):
    for label, replacement in MATH_LABELS.items():
        value = value.replace(label, replacement)
    return ax.text(x, y, value, fontsize=size, color=color,
                   fontproperties=FONT, va="top", **kwargs)


def eq(ax, x, y, value, size=27, color=BLUE):
    return text(ax, x, y, "$" + value + "$", size, color)


def page(title, subtitle, section):
    fig = plt.figure(figsize=(16, 9), facecolor=WHITE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    text(ax, .05, .955, section.upper(), 12, TEAL)
    text(ax, .05, .905, title, 31)
    text(ax, .05, .835, subtitle, 16, MUTED)
    ax.plot([.05, .95], [.79, .79], color=LINE, lw=1)
    number = len(PAGES) + 1
    text(ax, .05, .042, "XMVB-CPP  /  STRING → STRUCTURE  /  审阅初稿 v1 · 2026-10-07", 10, MUTED)
    text(ax, .95, .042, f"{number:02d}", 11, MUTED, ha="right")
    PAGES.append(fig)
    return ax


def box(ax, x, y, w, h, title, body="", color=BLUE):
    ax.add_patch(FancyBboxPatch((x, y-h), w, h,
                 boxstyle="round,pad=0.012,rounding_size=0.012",
                 ec=LINE, fc=PALE, lw=1))
    text(ax, x+.014, y-.016, title, 21, color)
    if body:
        text(ax, x+.014, y-.075, body, 17, linespacing=1.6)


def arrow(ax, start, end, color=MUTED):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                               mutation_scale=17, color=color, lw=1.8))


def takeaway(ax, value):
    ax.add_patch(Rectangle((.05, .087), .90, .09, fc="#EAF5F3", ec="none"))
    text(ax, .068, .152, value, 19, TEAL)


def table(ax, x, y, widths, rows, row_height=.078, size=17):
    for i, row in enumerate(rows):
        x0 = x
        for width, value in zip(widths, row):
            ax.add_patch(Rectangle((x0, y-(i+1)*row_height), width,
                         row_height, fc=INK if i == 0 else (PALE if i % 2 else WHITE),
                         ec=LINE, lw=.7))
            text(ax, x0+.012, y-i*row_height-.020, value, size,
                 WHITE if i == 0 else INK)
            x0 += width


def matrix(ax, x, y, values, row_labels=None, col_labels=None,
           cell=.074, font_size=20, highlight=True):
    rows, cols = len(values), len(values[0])
    if col_labels:
        for j, label in enumerate(col_labels):
            text(ax, x+(j+.5)*cell, y+.040, label, 15, MUTED, ha="center")
    for i in range(rows):
        if row_labels:
            text(ax, x-.014, y-(i+.30)*cell, row_labels[i], 15, MUTED, ha="right")
        for j in range(cols):
            v = values[i][j]
            active = highlight and str(v) not in ("0", "0.0", "")
            ax.add_patch(Rectangle((x+j*cell, y-(i+1)*cell), cell, cell,
                         fc="#DBECEF" if active else WHITE, ec=LINE, lw=1))
            text(ax, x+(j+.5)*cell, y-(i+.25)*cell, str(v), font_size,
                 TEAL if active else MUTED, ha="center")


def paired_terms(pairs):
    """Expand creation operators and account for canonical fermion parity.

    IDs 0..3 are alpha and 4..7 beta; both are ascending within spin.
    One common minus phase is chosen for the two teaching structures.
    """
    terms = {}
    for choices in itertools.product((0, 1), repeat=len(pairs)):
        operators, coefficient = [], -1
        for (i, j), swap in zip(pairs, choices):
            operators.extend((i, j+4) if swap == 0 else (i+4, j))
            coefficient *= 1 if swap == 0 else -1
        inversions = sum(operators[i] > operators[j]
                         for i in range(len(operators))
                         for j in range(i+1, len(operators)))
        coefficient *= (-1)**inversions
        ordered = sorted(operators)
        key = (tuple(v for v in ordered if v < 4),
               tuple(v-4 for v in ordered if v >= 4))
        terms[key] = terms.get(key, 0) + coefficient
    return terms


def verify():
    """Independently check expansion, direct contraction, and both adjoints."""
    left = paired_terms(((0, 1), (2, 3)))
    right = paired_terms(((0, 2), (1, 3)))
    expected_left = {((0, 2), (1, 3)): 1, ((0, 3), (1, 2)): 1,
                     ((1, 2), (0, 3)): 1, ((1, 3), (0, 2)): 1}
    expected_right = {((0, 1), (2, 3)): 1, ((0, 3), (1, 2)): -1,
                      ((1, 2), (0, 3)): -1, ((2, 3), (0, 1)): 1}
    assert left == expected_left and right == expected_right
    overlap = np.array([[1, .18, .07, .03], [.18, 1, .11, .06],
                        [.07, .11, 1, .16], [.03, .06, .16, 1]])
    assert np.linalg.eigvalsh(overlap).min() > 0
    alpha_l = [key[0] for key in left]
    beta_l = [key[1] for key in left]
    alpha_r = [key[0] for key in right]
    beta_r = [key[1] for key in right]
    sa = np.array([[np.linalg.det(overlap[np.ix_(i, j)])
                    for j in alpha_r] for i in alpha_l])
    sb = np.array([[np.linalg.det(overlap[np.ix_(i, j)])
                    for j in beta_r] for i in beta_l])
    cl = np.eye(4)
    cr = np.diag(list(right.values()))
    explicit = sum(left[kl]*right[kr] * sa[i, j]*sb[i, j]
                   for i, kl in enumerate(left) for j, kr in enumerate(right))
    regrouped = np.sum(cl*(sa@cr@sb.T))
    assert abs(explicit-regrouped) < 1e-13
    rng = np.random.default_rng(17)
    l, r = rng.normal(size=(3, 2)), rng.normal(size=(4, 5))
    a, b = rng.normal(size=(3, 4)), rng.normal(size=(2, 5))
    da, db = rng.normal(size=a.shape), rng.normal(size=b.shape)
    f = lambda aa, bb: float(np.sum(l*(aa@r@bb.T)))
    step = 1e-6
    grad_a, grad_b = l@b@r.T, l.T@a@r
    fd_a = (f(a+step*da, b)-f(a-step*da, b))/(2*step)
    fd_b = (f(a, b+step*db)-f(a, b-step*db))/(2*step)
    error_a = abs(fd_a-np.sum(grad_a*da))
    error_b = abs(fd_b-np.sum(grad_b*db))
    assert error_a < 1e-7 and error_b < 1e-7
    return {"scope": "document toy algebra only; not a VBSCF benchmark",
            "phase_convention": "common minus phase; sorted alpha then beta",
            "explicit_overlap": explicit, "regrouped_overlap": regrouped,
            "overlap_error": abs(explicit-regrouped),
            "adjoint_A_fd_error": error_a, "adjoint_B_fd_error": error_b,
            "cpp_revision": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "official_reference": "attachment xmvb4.0 branch a3613ed32ebffcaf2d98c1c9124555708bea995b"}


def render_content(checks):
    ax = page("从 spin-string 直接收缩到 VB structure",
              "主 H/S 收缩与 action：精确重排求和，跳过逐行列式对的组合与散射", "01 / 本次组会的主线")
    box(ax, .06, .72, .27, .27, "我们解决什么？",
        "结构展开产生大量项对\n完整 determinant-pair 的\n组合、查询和累加反复出现", OLD)
    box(ax, .365, .72, .27, .27, "我们怎样改变？",
        "共享 unique-spin 核\n直接进行矩阵／块收缩\n直接投影回所选结构空间")
    box(ax, .67, .72, .27, .27, "没有改变什么？",
        "固定轨道与积分核\n相同结构波函数\n收缩不引入额外近似", TEAL)
    text(ax, .07, .35, "先看一个两对 singlet 的例子，再讲公式与 Davidson。", 23)
    text(ax, .07, .27, "本稿只解释收缩算法；RI、正交 direct-CI、轨道优化不展开。", 18, MUTED)
    takeaway(ax, "方法贡献：改变数值收缩组织，而不是改写波函数。")

    ax = page("先区分三个层级", "展开项不是新的结构变分自由度", "02 / 表示")
    box(ax, .06, .71, .27, .24, "VB structure",
        "一个指定成键／自旋耦合\n的波函数构件\n最终求解结构系数 vᵢ", TEAL)
    box(ax, .365, .71, .27, .24, "完整 determinant",
        "一条 α 串与一条 β 串\n组合成的展开项\nD = (αₐ, βᵦ)", OLD)
    box(ax, .67, .71, .27, .24, "unique-spin-string",
        "去重后的单自旋占据串\n共享算子核的计算索引\n数量 Uα、Uβ", BLUE)
    arrow(ax, (.34, .59), (.35, .59))
    arrow(ax, (.645, .59), (.657, .59))
    eq(ax, .15, .38, r"|\Phi_I\rangle=\sum_{a,b}(C_I)_{ab}|a_\alpha b_\beta\rangle", 30)
    text(ax, .10, .28, "同一个 spin-string 可被许多展开项、许多结构共同使用。", 22)
    takeaway(ax, "保留必要的组合学关系，不等于保留 determinant-pair 数值计算层。")

    ax = page("两对 singlet：一个结构，四个展开项",
              "例子固定为 ΦA = (1,2)(3,4)；未归一化，统一选取整体负相位", "03 / 贯穿例子")
    eq(ax, .65, .73, r"|\Phi_A\rangle=-P_{12}^{\dagger}P_{34}^{\dagger}|0\rangle", 22)
    for center, labels in ((.23, (1, 2)), (.54, (3, 4))):
        ax.plot([center-.055, center+.055], [.69, .69], color=TEAL, lw=3)
        for x, label in ((center-.055, labels[0]), (center+.055, labels[1])):
            ax.scatter(x, .69, s=420, color=WHITE, edgecolor=TEAL, linewidth=2, zorder=3)
            text(ax, x, .713, str(label), 18, TEAL, ha="center")
    eq(ax, .10, .60, r"P_{ij}^{\dagger}=a_{i\alpha}^{\dagger}a_{j\beta}^{\dagger}-a_{i\beta}^{\dagger}a_{j\alpha}^{\dagger}", 24)
    table(ax, .10, .52, [.12, .22, .22, .24], [
        ["展开项", "α 占据串", "β 占据串", "系数（含重排符号）"],
        ["D₁", "(1,3)", "(2,4)", "+1"],
        ["D₂", "(1,4)", "(2,3)", "+1"],
        ["D₃", "(2,3)", "(1,4)", "+1"],
        ["D₄", "(2,4)", "(1,3)", "+1"],
    ], row_height=.063, size=16)
    text(ax, .10, .198, "约定：α 在前、β 在后，各自升序；选取一个共同整体负相位。", 13, MUTED)
    takeaway(ax, "两对 singlet → 4 个展开项；不能误画成只有 2 条 α / β strings。")

    ax = page("旧收缩：先组合完整行列式对，再累加结构",
              "XMVB4.0 已复用 string 核；剩下的瓶颈是 pair 组合和结构散射", "04 / 官方代码路径")
    matrix(ax, .125, .67, [[f"{i+1},{j+1}" for j in range(4)] for i in range(4)],
           ["D₁", "D₂", "D₃", "D₄"], ["D₁", "D₂", "D₃", "D₄"], cell=.066, font_size=16)
    text(ax, .075, .345, "单个结构自矩阵元：16 个有序项对\n实际代码还利用对称性与连接关系。", 16, MUTED, linespacing=1.7)
    for y, title, body in ((.71, "① Hamhd：单自旋核", "unique α / β 的 H、S"),
                            (.55, "② Hov1：组合 pair", "Ham = HαSβ + SαHβ + Hab"),
                            (.39, "③ Hovstr1：结构散射", "查归属、乘符号、累加 Hvb / Svb")):
        box(ax, .48, y, .44, .13, title, body, OLD if y != .71 else BLUE)
    text(ax, .06, .235, "官方没有保存完整 determinant H/S 矩阵；它仍形成 pair 标量并逐项累加。", 16)
    takeaway(ax, "我们消除的是数值收缩中间层，不是宣称首次发现 string 复用。")

    ax = page("新表示：把结构展开放进 string 系数矩阵",
              "行 = unique α strings；列 = unique β strings；非零元素 = 带符号的结构连接", "05 / 索引重排")
    text(ax, .08, .71, "左边的四个展开项", 22, OLD)
    for i, label in enumerate(("D₁ → α(1,3), β(2,4)", "D₂ → α(1,4), β(2,3)",
                               "D₃ → α(2,3), β(1,4)", "D₄ → α(2,4), β(1,3)")):
        text(ax, .085, .63-i*.071, label, 19)
    arrow(ax, (.43, .47), (.53, .47), TEAL)
    text(ax, .63, .72, "Cₐ：匹配列序下为对角矩阵", 21, TEAL)
    matrix(ax, .64, .62, np.eye(4, dtype=int).tolist(),
           ["α13", "α14", "α23", "α24"], ["β24", "β23", "β14", "β13"], cell=.063)
    text(ax, .08, .28, "一般结构保留带符号的稀疏连接；对角形状只是这个例子的特例。", 19)
    text(ax, .08, .218, "不同结构局部列表可不同；核的行、列必须按左右 support 对齐。", 17, MUTED)
    takeaway(ax, "Cᵢ 不是新的近似波函数，而是原结构展开的精确重标记。")

    ax = page("重叠：同一个求和，换一种计算顺序",
              "先吸收一侧自旋，再吸收另一侧；最后直接投影到结构", "06 / 核心公式")
    eq(ax, .07, .72, r"S_{IJ}=\sum_{a,a',b,b'}(C_I)_{ab}S^\alpha_{aa'}(C_J)_{a'b'}S^\beta_{bb'}", 26)
    arrow(ax, (.50, .62), (.50, .56), TEAL)
    eq(ax, .20, .53, r"S_{IJ}=\langle C_I,\ S^\alpha C_J(S^\beta)^T\rangle_F", 32)
    table(ax, .09, .41, [.12, .40, .30], [
        ["步骤", "计算", "结果形状（局部 support）"],
        ["1", "W = Cⱼ (Sβ)ᵀ", "αⱼ × βᵢ"],
        ["2", "Y = Sα W", "αᵢ × βᵢ"],
        ["3", "Sᵢⱼ = Σₐᵦ (Cᵢ)ₐᵦ Yₐᵦ", "一个结构矩阵元"],
    ], row_height=.053, size=16)
    text(ax, .10, .185, "实数约定：〈L,R〉F = Tr(LᵀR)。例子中的稀疏 C 还可用专门收缩。", 13, MUTED)
    takeaway(ax, "中间只有 string-space 图像，不再形成并散射每个完整 determinant-pair 标量。")

    ax = page("Hamiltonian：沿用同一个收缩框架",
              "算子核可以不同，但两端始终是结构系数块 Cᵢ 与 Cⱼ", "07 / 多个物理通道")
    text(ax, .07, .72, "一电子 + same-spin 双电子", 22, TEAL)
    eq(ax, .09, .655, r"H^{\rm same}_{IJ}=\langle C_I,\ H^\alpha C_J(S^\beta)^T+S^\alpha C_J(H^\beta)^T\rangle_F", 24)
    text(ax, .09, .575, "Hα、Hβ 含有效一电子和各自自旋内的双电子贡献。", 17, MUTED)
    text(ax, .07, .50, "Opposite-spin 双电子", 22, TEAL)
    eq(ax, .09, .435, r"H^{\rm opp}_{IJ}=\sum_{p,q}V_{pq}\langle C_I,\ D^{\alpha,p}C_J(D^{\beta,q})^T\rangle_F", 25)
    text(ax, .09, .35, "p、q 是轨道对通道；D 是未除以 overlap 的单自旋跃迁密度核。\nExact / RI 只改变积分核的表示与应用方式，不改变结构投影。", 17, MUTED, linespacing=1.7)
    text(ax, .09, .235, "不活性空间的标量项另加 Ecore Sᵢⱼ；本页不省略其物理作用。", 16)
    takeaway(ax, "不仅 overlap：各 Hamiltonian 通道都可直接收缩到 structure。")

    ax = page("Backward：反向权重也在 string 空间收缩",
              "从一条通用双线性收缩看梯度，再对它做方向微分得到 HVP 所需响应", "08 / 导数闭合")
    eq(ax, .20, .72, r"F=\langle L,\ A R B^T\rangle_F", 32)
    box(ax, .075, .59, .38, .20, "给 α 核的反向权重", "先把 β 核吸收到系数两端", TEAL)
    box(ax, .545, .59, .38, .20, "给 β 核的反向权重", "先把 α 核吸收到系数两端", BLUE)
    eq(ax, .15, .445, r"\frac{\partial F}{\partial A}=LBR^T", 27, TEAL)
    eq(ax, .61, .445, r"\frac{\partial F}{\partial B}=L^TAR", 27)
    text(ax, .085, .31, "不必先求每个完整 determinant pair 的导数，再散射到结构。", 21)
    text(ax, .085, .245, "VBSCF 能量梯度还要组合 dH − E dS；HVP 还含核及结构系数的响应。", 16, MUTED)
    takeaway(ax, "这两条是收缩 adjoint，不是完整轨道梯度／relaxed Hessian 的全部公式。")

    ax = page("Davidson：搜索空间与计算载体不是一回事",
              "我们的未知量是结构系数；string tensor 只承担内部 action", "09 / 求解空间")
    eq(ax, .07, .72, r"H_{\rm str}v=E S_{\rm str}v,\qquad v\in\mathbb{R}^{N_{\rm str}}", 28)
    for x, w, title, body in ((.065, .24, "结构向量 v", "展开 X = Σᵢ vᵢ Cᵢ"),
                               (.38, .24, "string action", "得到 YH、YS"),
                               (.695, .24, "结构 image", "投影 (Hv)ᵢ = 〈Cᵢ,YH〉")):
        box(ax, x, .59, w, .15, title, body, TEAL)
    arrow(ax, (.32, .51), (.36, .51))
    arrow(ax, (.635, .51), (.675, .51))
    table(ax, .075, .37, [.32, .26, .27], [
        ["已核对的实际模块", "迭代／求解空间", "求解方法"],
        ["XMVB4.0 普通 RDM–VBSCF", "VB structure", "dense 广义对角化"],
        ["XMVB4.0 的 FCI 模块", "完整 determinant", "Davidson"],
        ["我们的 VBSCF Davidson", "所选 VB structure", "matrix-free Davidson"],
    ], row_height=.050, size=15)
    takeaway(ax, "FCI 的 determinant Davidson 没有错；不能把它与 selected-VB 问题混为一谈。")

    ax = page("贡献在哪里，边界又在哪里？",
              "先建立准确的方法主张，再用隔离 benchmark 衡量它的收益", "10 / 总结与审阅")
    table(ax, .06, .72, [.23, .34, .31], [
        ["问题", "我们已经做到", "不能据此声称"],
        ["主 H/S 收缩与 action", "跳过逐 determinant-pair 散射", "所有辅助路径都跳过 pair"],
        ["核与导数复用", "矩阵／块收缩及 string adjoint", "任意 pair 只计算一次"],
        ["Davidson", "结构空间迭代，不组装完整 H/S", "内部 carrier 没有大维度"],
        ["Tile", "限制 pair 载荷与局部 workspace", "整个进程内存严格有界"],
    ], row_height=.065, size=15)
    eq(ax, .10, .355, r"\mathrm{cost/channel/vector}\sim U_\alpha^2U_\beta+U_\alpha U_\beta^2", 26)
    text(ax, .08, .284, "完整稠密乘积、单可分离通道的 action：U⁴ 主循环 → U³。\n不计核生成、通道数与整体梯度/HVP 成本；U 本身仍可能指数增长。\n单结构的稀疏特例不能直接套这个加速比。", 16, MUTED, linespacing=1.45)
    takeaway(ax, "重要突破是收缩组织；不是已经证明 VBSCF 对活性空间大小变成多项式算法。")

    ax = page("附录 A：用第二个结构检查公式，不只看自重叠",
              "ΦB = (1,3)(2,4)；同一整体相位约定，系数为 +1、−1、−1、+1", "附录 / 可复现教学检查")
    table(ax, .07, .72, [.18, .20, .20, .28], [
        ["结构 B 展开项", "α 占据串", "β 占据串", "系数"],
        ["1", "(1,2)", "(3,4)", "+1"],
        ["2", "(1,4)", "(2,3)", "−1"],
        ["3", "(2,3)", "(1,4)", "−1"],
        ["4", "(3,4)", "(1,2)", "+1"],
    ], row_height=.058, size=16)
    text(ax, .08, .39, "取一个固定、正定、非正交的 4×4 单电子 overlap，\n分别用 16 项直接求和与矩阵收缩计算未归一化的 Sₐᵦ。", 17, linespacing=1.65)
    table(ax, .08, .29, [.36, .48], [
        ["检查", "脚本运行结果"],
        ["直接和 / 重排值", f"{checks['explicit_overlap']:.10f} / {checks['regrouped_overlap']:.10f}"],
    ], row_height=.047, size=16)
    takeaway(ax, f"代数一致性误差 {checks['overlap_error']:.2e}；仅为教学检查，不是实际分子测速。")

    ax = page("附录 B：代码证据与下一轮审阅问题",
              "来源固定到当前 C++ 快照与附件 xmvb4.0 分支；不使用 DeepVBH", "附录 / 证据索引")
    table(ax, .06, .72, [.31, .57], [
        ["命题", "来源文件与行号"],
        ["官方 pair → structure 散射", "src/vbscf/gradient_rdm.F90:3380  Hovstr1"],
        ["官方普通 VBSCF dense 求解", "src/math/eigencalc.c:85  LAPACKE_dsygv"],
        ["官方 FCI Davidson = det_num", "src/mc_scf/fci_cal.c:467  init_davidson"],
        ["我们的结构矩阵收缩", "src/vbscf/structures/assembly/block_kernels.hpp:31"],
        ["我们的结构 expand / project", "src/vbscf/structures/assembly/action.cpp:858 / 778"],
        ["当前仍有初始化展开", "src/vbscf/structures/expansion/expander.cpp:66"],
    ], row_height=.050, size=14)
    text(ax, .07, .31, "下一轮：例子的层级是否直观？符号约定是否清楚？\n贡献表述是否准确？应补哪组同算例、同精度的隔离测速？", 18, linespacing=1.7)
    text(ax, .07, .21, "完整版本号与教学检查数值见同目录 string_contraction_checks.json。", 14, MUTED)
    takeaway(ax, "先共同审阅第 3–6 页的例子与主公式，再逐轮加深和补充证据。")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"docs/string_contraction.pdf")
    parser.add_argument("--font", type=Path,
                        default=Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"))
    args = parser.parse_args()
    global FONT
    if not args.font.is_file():
        raise SystemExit(f"CJK font not found: {args.font}; specify --font")
    FONT = FontProperties(fname=str(args.font))
    # Type 3 embeds vector glyphs without mislabelling CFF outlines as TrueType.
    plt.rcParams.update({"pdf.fonttype": 3, "ps.fonttype": 3,
                         "mathtext.fontset": "stix", "axes.unicode_minus": False})
    checks = verify()
    checks["render_versions"] = {"matplotlib": matplotlib.__version__,
                                 "numpy": np.__version__, "font": str(args.font)}
    render_content(checks)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    thumbnails = []
    with PdfPages(args.output, metadata={
        "Title": "VBSCF: exact spin-string-to-structure contraction — review draft v1",
        "Author": "xmvb-cpp project", "Subject": "Chinese group-meeting teaching deck"}) as pdf:
        for fig in PAGES:
            pdf.savefig(fig)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=45, facecolor=WHITE)
            buf.seek(0)
            thumbnails.append(Image.open(buf).convert("RGB"))
            plt.close(fig)
    width, height = thumbnails[0].size
    contact = Image.new("RGB", (3*(width+12)+12, 4*(height+12)+12), "#D7E1EB")
    for i, thumb in enumerate(thumbnails):
        contact.paste(ImageOps.expand(thumb, border=1, fill=LINE),
                      (12+(i % 3)*(width+12), 12+(i // 3)*(height+12)))
    preview = args.output.with_name(args.output.stem+"_preview.png")
    contact.save(preview)
    checks["pages"] = len(PAGES)
    checks_path = args.output.with_name(args.output.stem+"_checks.json")
    checks_path.write_text(json.dumps(checks, ensure_ascii=False, indent=2)+"\n")
    print(args.output)
    print(preview)
    print(json.dumps(checks, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

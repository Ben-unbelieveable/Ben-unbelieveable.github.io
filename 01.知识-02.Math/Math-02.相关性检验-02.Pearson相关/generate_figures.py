#!/usr/bin/env python3
"""生成 Pearson 相关文档配图。"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib import font_manager
from scipy import stats

OUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(42)

# 尝试加载常见中文字体
for name in ["PingFang SC", "Heiti SC", "STHeiti", "Arial Unicode MS", "SimHei"]:
    if any(name in f.name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
        break
plt.rcParams["axes.unicode_minus"] = False


def save(fig, name: str) -> None:
    """保存图片到当前目录。"""
    fig.savefig(OUT / name, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig01_intuition() -> None:
    """图1：不同 |r| 对应的散点形态。"""
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.2))
    specs = [
        ("强正相关", 0.92, 1.0),
        ("弱相关", 0.35, 1.2),
        ("无相关", 0.0, 1.0),
        ("强负相关", -0.88, 1.0),
    ]
    for ax, (title, rho, noise) in zip(axes, specs):
        x = RNG.normal(0, 1, 60)
        y = rho * x + RNG.normal(0, noise, 60)
        r = np.corrcoef(x, y)[0, 1]
        ax.scatter(x, y, s=22, alpha=0.75, c="#2563eb", edgecolors="white", linewidths=0.4)
        if abs(rho) > 0.5:
            m, b = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 50)
            ax.plot(xs, m * xs + b, color="#dc2626", lw=1.8, alpha=0.9)
        ax.axhline(0, color="#e5e7eb", lw=0.8)
        ax.axvline(0, color="#e5e7eb", lw=0.8)
        ax.set_title(f"{title}\n$r$ = {r:.2f}", fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#d1d5db")
    fig.suptitle("Pearson $r$：散点越接近直线，$|r|$ 越大", fontsize=13, y=1.02)
    fig.text(0.5, -0.02, "圆形云团 → $r \\approx 0$（无线性趋势）", ha="center", fontsize=10, color="#4b5563")
    save(fig, "fig01-intuition.png")


def fig02_scatter_and_test() -> None:
    """图2：公式直观 + t 检验示意。"""
    fig = plt.figure(figsize=(11, 4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])

    ax1 = fig.add_subplot(gs[0, 0])
    n = 50
    x = RNG.normal(0, 1, n)
    y = 0.8 * x + RNG.normal(0, 0.55, n)
    r, p = stats.pearsonr(x, y)
    ax1.scatter(x, y, s=30, c="#0d9488", alpha=0.8, edgecolors="white", linewidths=0.5)
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ax1.plot(xs, m * xs + b, color="#b45309", lw=2, label="最小二乘拟合")
    ax1.set_xlabel("$X$")
    ax1.set_ylabel("$Y$")
    ax1.set_title(f"样本相关 $r$ = {r:.3f}", fontsize=12)
    ax1.legend(fontsize=9, loc="upper left")

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.axis("off")
    formula = (
        r"$r = \dfrac{\sum (x_i-\bar{x})(y_i-\bar{y})}"
        r"{\sqrt{\sum(x_i-\bar{x})^2}\,\sqrt{\sum(y_i-\bar{y})^2}}$"
        "\n\n"
        r"$H_0\!: \rho = 0 \quad\Rightarrow\quad "
        r"t = \dfrac{r\sqrt{n-2}}{\sqrt{1-r^2}} \sim t_{n-2}$"
        f"\n\n本例：$n={n}$，$t \\approx {r * np.sqrt(n-2) / np.sqrt(1-r**2):.2f}$，$p = {p:.2e}$"
    )
    ax2.text(
        0.05, 0.95, formula, va="top", fontsize=12,
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#f8fafc", edgecolor="#cbd5e1"),
    )
    ax2.text(
        0.05, 0.12,
        "分子：$X,Y$ 同向偏离均值的乘积之和\n"
        "分母：两变量各自离差的标准化尺度",
        va="bottom", fontsize=10, color="#374151",
    )
    fig.suptitle("Pearson 系数与显著性检验", fontsize=13, y=1.02)
    save(fig, "fig02-pmf.png")


def fig03_scenarios() -> None:
    """图3：典型适用场景示意。"""
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))

    # 酶活 vs 荧光
    ax = axes[0]
    enzyme = RNG.uniform(10, 100, 40)
    fluo = 2.1 * enzyme + RNG.normal(0, 15, 40)
    r1 = np.corrcoef(enzyme, fluo)[0, 1]
    ax.scatter(enzyme, fluo, s=28, c="#7c3aed", alpha=0.8)
    ax.set_xlabel("酶活 (U/mL)")
    ax.set_ylabel("荧光强度 (a.u.)")
    ax.set_title(f"酶活 vs 荧光\n$r$ = {r1:.2f}")

    # 两基因表达
    ax = axes[1]
    g1 = RNG.lognormal(2, 0.6, 45)
    g2 = 0.7 * g1 + RNG.lognormal(2, 0.5, 45) * 0.4
    r2 = np.corrcoef(np.log1p(g1), np.log1p(g2))[0, 1]
    ax.scatter(np.log1p(g1), np.log1p(g2), s=28, c="#059669", alpha=0.8)
    ax.set_xlabel("基因 A (log1p)")
    ax.set_ylabel("基因 B (log1p)")
    ax.set_title(f"共表达基因\n$r$ = {r2:.2f}")

    # biomarker 初筛
    ax = axes[2]
    bio = RNG.normal(0, 1, 50)
    pheno = 0.65 * bio + RNG.normal(0, 0.9, 50)
    r3 = np.corrcoef(bio, pheno)[0, 1]
    ax.scatter(bio, pheno, s=28, c="#ea580c", alpha=0.8)
    ax.set_xlabel("Biomarker 水平")
    ax.set_ylabel("表型评分")
    ax.set_title(f"Biomarker 初筛\n$r$ = {r3:.2f}，$R^2$ = {r3**2:.2f}")

    fig.suptitle("Pearson 适用：连续变量、近似线性关系", fontsize=13, y=1.02)
    save(fig, "fig03-scenarios.png")


def fig04_limitations() -> None:
    """图4：Anscombe 型局限示意。"""
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7))
    panels = []

    # 正常线性
    x = RNG.normal(5, 1.5, 20)
    y = 0.5 * x + 3 + RNG.normal(0, 0.8, 20)
    panels.append((x, y, "线性（适用）"))

    # 离群点
    x2 = x.copy()
    y2 = y.copy()
    x2[-1], y2[-1] = 12, 2
    panels.append((x2, y2, "离群点拉高 $r$"))

    # U 形
    x3 = RNG.uniform(-2, 2, 30)
    y3 = x3 ** 2 + RNG.normal(0, 0.3, 30)
    panels.append((x3, y3, "U 形：$r \\approx 0$"))

    # 强相关但离群
    x4 = RNG.normal(5, 1.2, 20)
    y4 = 0.5 * x4 + 3 + RNG.normal(0, 0.6, 20)
    x4[-1], y4[-1] = 9, 2
    panels.append((x4, y4, "同一 $r$，形态不同"))

    for ax, (px, py, title) in zip(axes.flat, panels):
        r = np.corrcoef(px, py)[0, 1]
        color = "#dc2626" if "离群" in title or "U" in title else "#2563eb"
        ax.scatter(px, py, s=32, c=color, alpha=0.85, edgecolors="white", linewidths=0.4)
        ax.set_title(f"{title}\n$r$ = {r:.2f}", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#d1d5db")

    fig.suptitle("局限：必画散点图；非线性 / 离群 → 考虑 Spearman 或变换", fontsize=12, y=1.01)
    fig.text(0.5, 0.01, "相关 $\\neq$ 因果：需实验设计验证", ha="center", fontsize=10, color="#6b7280")
    save(fig, "fig04-limitations.png")


if __name__ == "__main__":
    fig01_intuition()
    fig02_scatter_and_test()
    fig03_scenarios()
    fig04_limitations()
    print("done:", sorted(p.name for p in OUT.glob("fig*.png")))

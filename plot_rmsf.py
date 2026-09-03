# -*- coding: utf-8 -*-
"""
plot_rmsf.py — 绘制 VMD RMSF 计算结果曲线图
=================================================
读取 rmsf_analysis.tcl 生成的 .residue.dat 或 .atom.dat 文件，
绘制 残基编号-RMSF 曲线，高亮柔性区域（> mean+1σ），输出 PNG/PDF。

用法:
    python plot_rmsf.py rmsf_result.residue.dat
    python plot_rmsf.py rmsf_result.residue.dat -o myfig --title "Protein X"
    python plot_rmsf.py rmsf_result.residue.dat --sigma 1.5 --no-highlight

依赖: matplotlib (pip install matplotlib)
"""
import argparse
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")  # 无显示环境也可出图
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
except ImportError:
    sys.exit("缺少 matplotlib，请先安装: pip install matplotlib")


def setup_chinese_font() -> str:
    """配置中文字体（Windows 优先微软雅黑/黑体，失败则回退英文标签）。"""
    candidates = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "PingFang SC", "Arial"]
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in available:
            plt.rcParams["font.sans-serif"] = [name]
            plt.rcParams["axes.unicode_minus"] = False
            return name
    return ""


def read_dat(path: str):
    """解析 dat 文件 -> (resid 列表, rmsf 列表, header 元信息字符串)。"""
    resids, rmsfs = [], []
    header = []
    ncols = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header.append(line.lstrip("# ").strip())
                continue
            parts = line.split()
            if ncols is None:
                ncols = len(parts)
            # 两种格式: residue.dat(5列: resid resname chain segid rmsf)
            #           atom.dat   (6列: ... atomname rmsf)
            if ncols == 5:
                rid, rmsf = parts[0], parts[4]
            elif ncols == 6:
                rid, rmsf = parts[0], parts[5]
            else:
                sys.exit(f"无法识别的数据列数 ({ncols}): {line}")
            try:
                resids.append(int(rid))
                rmsfs.append(float(rmsf))
            except ValueError:
                continue
    if not resids:
        sys.exit(f"文件中未解析到有效数据: {path}")
    return resids, rmsfs, " | ".join(header)


def main():
    ap = argparse.ArgumentParser(description="绘制 VMD RMSF 结果曲线")
    ap.add_argument("dat", help="rmsf_analysis.tcl 输出的 .residue.dat 或 .atom.dat")
    ap.add_argument("-o", "--output", default=None, help="输出文件前缀（默认与输入同名）")
    ap.add_argument("--title", default=None, help="图标题")
    ap.add_argument("--sigma", type=float, default=1.0,
                    help="柔性区域阈值 = mean + N*sigma（默认 1.0）")
    ap.add_argument("--no-highlight", action="store_true", help="不高亮柔性区域")
    args = ap.parse_args()

    font = setup_chinese_font()
    resids, rmsfs, header = read_dat(args.dat)

    n = len(resids)
    mean = sum(rmsfs) / n
    std = (sum((x - mean) ** 2 for x in rmsfs) / n) ** 0.5
    thresh = mean + args.sigma * std
    flexible = [i for i, v in enumerate(rmsfs) if v > thresh]

    prefix = args.output or str(Path(args.dat).with_suffix(""))
    title = args.title or (header.split("|")[0].strip() if header else "RMSF")

    # ---------------- 绘图 ----------------
    fig, ax = plt.subplots(figsize=(12, 5), dpi=150)
    ax.plot(resids, rmsfs, color="#1565C0", lw=1.3, marker="o", ms=2.5,
            mfc="#1565C0", mec="none", label="RMSF", zorder=3)

    if not args.no_highlight and flexible:
        fx = [resids[i] for i in flexible]
        fy = [rmsfs[i] for i in flexible]
        ax.scatter(fx, fy, color="#D32F2F", s=22, zorder=4,
                   label=f"柔性残基 (>{mean:.2f}+{args.sigma:g}σ={thresh:.2f} Å, n={len(flexible)})")

    ax.axhline(mean, color="#2E7D32", ls="--", lw=1.0,
               label=f"均值 = {mean:.2f} Å")
    ax.axhline(thresh, color="#EF6C00", ls=":", lw=1.0,
               label=f"阈值 = {thresh:.2f} Å")

    ax.set_xlabel("残基编号 Residue ID" if font else "Residue ID")
    ax.set_ylabel("RMSF (Å)" if font else "RMSF ($\\AA$)")
    ax.set_title(f"{title} — RMSF per Residue" if font else f"{title} — RMSF per Residue")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.3, ls="--")
    fig.tight_layout()

    png, pdf = f"{prefix}.png", f"{prefix}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    # ---------------- 摘要 ----------------
    imax = max(range(n), key=lambda i: rmsfs[i])
    print(f"[plot_rmsf] 读取 {n} 个数据点，来自: {args.dat}")
    print(f"[plot_rmsf] 均值={mean:.3f} Å, 标准差={std:.3f} Å, "
          f"阈值(mean+{args.sigma:g}σ)={thresh:.3f} Å")
    print(f"[plot_rmsf] 柔性残基数: {len(flexible)} / {n}"
          + (f"，例如: {[resids[i] for i in flexible[:10]]}" if flexible else ""))
    print(f"[plot_rmsf] 最大RMSF残基: {resids[imax]} ({rmsfs[imax]:.3f} Å)")
    print(f"[plot_rmsf] 已输出: {png} , {pdf}")


if __name__ == "__main__":
    main()

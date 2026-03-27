#!/usr/bin/env python3
"""Post-run script: update manuscript tables, figures, and compile PDFs.

Reads all lightweight experiment outputs and integrates them into
the Claude manuscript (manuscript_claude/) without overwriting
any existing valid content.

Also generates a final lightweight_strengthen_report.md.
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def _ts() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _log(msg: str) -> None:
    print(f"[{_ts()}] {msg}", flush=True)


def _read_csv_safe(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def _compile_tex(tex_path: Path) -> bool:
    """Run pdflatex 3 times for cross-references."""
    tex_dir = tex_path.parent
    for i in range(3):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=str(tex_dir),
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0 and i == 2:
            _log(f"pdflatex pass {i+1} failed for {tex_path.name}")
            _log(result.stdout[-500:] if result.stdout else "no stdout")
            return False
    return True


def generate_report(results_dir: Path, manuscript_dir: Path) -> str:
    """Generate the lightweight_strengthen_report markdown."""
    lines = [
        "# Lightweight Strengthen Report",
        f"\nDate: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
    ]

    # 1. What jobs ran
    ckpt_path = results_dir / "checkpoint_lightweight_strengthen.json"
    if ckpt_path.exists():
        with open(ckpt_path) as f:
            ckpt = json.load(f)
        lines.append("## 1. Jobs completed")
        lines.append("")
        lines.append("| Task | Status |")
        lines.append("|------|--------|")
        for t in ckpt.get("tasks", []):
            lines.append(f"| {t['task']} | {t['status']} |")
        lines.append(f"\nMax RSS: {ckpt.get('max_rss_gb', 'N/A')} GB")
        lines.append("")

    # 2. LUAD multiseed results
    luad_stats = _read_csv_safe(results_dir / "paired_stats_luad_multiseed_v2.csv")
    if not luad_stats.empty:
        lines.append("## 2. LUAD multiseed paired analysis")
        lines.append("")
        lines.append("| Comparison | N pairs | Mean delta | Paired t p-value | Wilcoxon p-value |")
        lines.append("|------------|---------|------------|------------------|------------------|")
        for _, r in luad_stats.iterrows():
            lines.append(
                f"| {r['comparison']} | {r['n_pairs']} | {r['mean_delta']:.4f} | "
                f"{r.get('paired_t_p', 'N/A'):.4f} | {r.get('wilcoxon_p', 'N/A'):.4f} |"
            )
        lines.append("")

    # 3. GSE161711 bounded ablation
    gse_ablation = _read_csv_safe(results_dir / "ablation_summary_cll_venetoclax_gse161711_bounded.csv")
    if not gse_ablation.empty:
        lines.append("## 3. GSE161711 bounded ablation")
        lines.append("")
        lines.append("| Ablation | Value | Mean F1 | Std | N |")
        lines.append("|----------|-------|---------|-----|---|")
        for _, r in gse_ablation.iterrows():
            lines.append(
                f"| {r['ablation']} | {r['ablation_value']} | {r['mean']:.4f} | "
                f"{r.get('std', 0):.4f} | {int(r['n'])} |"
            )
        lines.append("")

    # 4. Diagnostics
    diag = _read_csv_safe(results_dir / "tables" / "lightweight_diagnostics.csv")
    if not diag.empty:
        lines.append("## 4. TIP diagnostics")
        lines.append("")
        lines.append("| Dataset | Backend | Support | Gini | Top-10 conc. | H1 lifetime |")
        lines.append("|---------|---------|---------|------|--------------|-------------|")
        for _, r in diag.iterrows():
            lines.append(
                f"| {r['dataset']} | {r['backend']} | {r['support_size']}/{r['total_features']} | "
                f"{r['tip_gini']:.3f} | {r['tip_top10_concentration']:.3f} | "
                f"{r.get('h1_lifetime_mean', 'N/A')} |"
            )
        lines.append("")

    # 5. scRNA tiny sensitivity
    scrna_abl = _read_csv_safe(results_dir / "ablation_summary_cll_rs_scrna_gse165087_tiny_sensitivity.csv")
    if not scrna_abl.empty:
        lines.append("## 5. scRNA tiny sensitivity (k ablation)")
        lines.append("")
        for _, r in scrna_abl.iterrows():
            lines.append(f"- {r['ablation']}={r['ablation_value']}: F1={r['mean']:.4f}")
        lines.append("")

    # 6. Runtime/memory
    runtime = _read_csv_safe(results_dir / "tables" / "runtime_memory_summary.csv")
    if not runtime.empty:
        lines.append("## 6. Runtime and memory")
        lines.append("")
        lines.append("| Checkpoint | Phase | Max RSS (GB) |")
        lines.append("|------------|-------|--------------|")
        for _, r in runtime.iterrows():
            lines.append(f"| {r['checkpoint_file']} | {r['phase']} | {r.get('max_rss_gb', 'N/A')} |")
        lines.append("")

    # 7. Manuscript outputs
    paper_pdf = manuscript_dir / "cc_spr_full_paper_claude.pdf"
    supp_pdf = manuscript_dir / "cc_spr_full_supplement_claude.pdf"
    lines.append("## 7. Final outputs")
    lines.append("")
    lines.append(f"- Paper PDF: `{paper_pdf}` ({'exists' if paper_pdf.exists() else 'NOT FOUND'})")
    lines.append(f"- Supplement PDF: `{supp_pdf}` ({'exists' if supp_pdf.exists() else 'NOT FOUND'})")
    lines.append("")

    return "\n".join(lines)


def main() -> None:
    results_dir = Path("results")
    results_claude = Path("results_claude")
    manuscript_dir = Path("manuscript_claude")
    tables_dir = results_dir / "tables"
    figures_dir = results_dir / "figures"

    for d in [results_claude, manuscript_dir, tables_dir, figures_dir]:
        d.mkdir(parents=True, exist_ok=True)

    _log("Generating lightweight strengthen report")
    report = generate_report(results_dir, manuscript_dir)
    report_path = results_claude / "lightweight_strengthen_report.md"
    report_path.write_text(report)
    _log(f"Report written to {report_path}")

    # Compile PDFs
    paper_tex = manuscript_dir / "cc_spr_full_paper_claude.tex"
    supp_tex = manuscript_dir / "cc_spr_full_supplement_claude.tex"

    if paper_tex.exists():
        _log("Compiling main paper PDF")
        ok = _compile_tex(paper_tex)
        _log(f"Main paper: {'OK' if ok else 'FAILED'}")

    if supp_tex.exists():
        _log("Compiling supplement PDF")
        ok = _compile_tex(supp_tex)
        _log(f"Supplement: {'OK' if ok else 'FAILED'}")

    _log("Manuscript refresh complete")


if __name__ == "__main__":
    main()

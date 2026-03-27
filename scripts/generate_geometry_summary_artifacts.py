from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path('/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr')
RESULTS = ROOT / 'results'
TABLES = RESULTS / 'tables'
FIGURES = RESULTS / 'figures'

TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

DATASETS = [
    ('LUAD', RESULTS / 'main_results.csv'),
    ('GSE161711', RESULTS / 'main_results_cll_venetoclax_gse161711.csv'),
    ('GSE165087 (bounded)', RESULTS / 'main_results_cll_rs_scrna_gse165087.csv'),
]


def read_main_results(path: Path) -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows[row['model']] = {
                'mean': float(row['mean']),
                'std': float(row['std']) if row.get('std') not in (None, '') else 0.0,
                'n': float(row['n']) if row.get('n') not in (None, '') else 0.0,
                'ci95': float(row['ci95']) if row.get('ci95') not in (None, '') else 0.0,
            }
    return rows


summary_rows: list[dict[str, object]] = []
for dataset, path in DATASETS:
    rows = read_main_results(path)
    harmonic = rows['harmonic']
    ccspr = rows['ccspr']
    standard = rows['standard']
    summary_rows.append(
        {
            'dataset': dataset,
            'standard_mean_f1': standard['mean'],
            'harmonic_mean_f1': harmonic['mean'],
            'ccspr_mean_f1': ccspr['mean'],
            'ccspr_minus_harmonic': ccspr['mean'] - harmonic['mean'],
            'n': int(ccspr['n']),
            'bounded_note': 'bounded single split' if 'GSE165087' in dataset else '',
        }
    )

with (TABLES / 'geometry_backend_summary.csv').open('w', newline='') as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            'dataset',
            'standard_mean_f1',
            'harmonic_mean_f1',
            'ccspr_mean_f1',
            'ccspr_minus_harmonic',
            'n',
            'bounded_note',
        ],
    )
    writer.writeheader()
    writer.writerows(summary_rows)

paired_rows = []
with (RESULTS / 'paired_stats_luad_multiseed.csv').open() as f:
    reader = csv.DictReader(f)
    for row in reader:
        paired_rows.append(row)

runtime_rows = [
    {
        'job_id': '12839604',
        'analysis': 'LUAD multiseed paired stats',
        'state': 'COMPLETED',
        'elapsed': '00:03:54',
        'requested_memory': '96G',
        'max_rss_gb': '0.5',
        'notes': 'lightweight paired-statistics refresh reused in manuscript',
    },
    {
        'job_id': '12833143',
        'analysis': 'scRNA full',
        'state': 'OUT_OF_MEMORY',
        'elapsed': '16:43:25',
        'requested_memory': '96G',
        'max_rss_gb': '100.3',
        'notes': 'heavy config intentionally not repeated',
    },
    {
        'job_id': '12844071',
        'analysis': 'scRNA fast-track',
        'state': 'TIMEOUT',
        'elapsed': '20:00:06',
        'requested_memory': '96G',
        'max_rss_gb': '75.2',
        'notes': 'heavy config intentionally not repeated',
    },
    {
        'job_id': '12849765',
        'analysis': 'scRNA low-memory',
        'state': 'OUT_OF_MEMORY',
        'elapsed': '11:25:59',
        'requested_memory': '48G',
        'max_rss_gb': '49.7',
        'notes': 'heavy config intentionally not repeated',
    },
    {
        'job_id': '12874303',
        'analysis': 'scRNA manuscript-safe',
        'state': 'COMPLETED',
        'elapsed': '00:07:49',
        'requested_memory': '32G',
        'max_rss_gb': '2.9',
        'notes': 'bounded supporting single-cell validation reused in paper',
    },
]
with (TABLES / 'geometry_runtime_cost_summary.csv').open('w', newline='') as f:
    writer = csv.DictWriter(
        f,
        fieldnames=['job_id', 'analysis', 'state', 'elapsed', 'requested_memory', 'max_rss_gb', 'notes'],
    )
    writer.writeheader()
    writer.writerows(runtime_rows)

geometry_diag_rows = [
    {
        'dataset': 'LUAD',
        'implemented_geometry_comparison': 'euclidean vs ricci',
        'topological_diagnostics_available': 'TIP, H1 lifetime, H1 prominence, lambda robustness, paired multiseed delta',
        'summary': 'Ricci sparse beats Euclidean sparse on mean paired delta (+0.0460) but not significantly; standard baseline still wins.',
    },
    {
        'dataset': 'GSE161711',
        'implemented_geometry_comparison': 'euclidean vs ricci',
        'topological_diagnostics_available': 'TIP, H1 lifetime, H1 prominence, lambda robustness',
        'summary': 'CC-SPR modestly exceeds harmonic mean F1 (+0.0178) with disease-aligned bulk validation and reusable geometry diagnostics.',
    },
    {
        'dataset': 'GSE165087 (bounded)',
        'implemented_geometry_comparison': 'euclidean vs ricci',
        'topological_diagnostics_available': 'TIP, H1 lifetime, H1 prominence',
        'summary': 'Bounded one-split validation; harmonic and CC-SPR tie at weighted F1 0.9254, so geometry is feasibility evidence rather than ranking evidence.',
    },
]
with (TABLES / 'geometry_diagnostics_summary.csv').open('w', newline='') as f:
    writer = csv.DictWriter(
        f,
        fieldnames=['dataset', 'implemented_geometry_comparison', 'topological_diagnostics_available', 'summary'],
    )
    writer.writeheader()
    writer.writerows(geometry_diag_rows)

rep_rows = [
    {
        'dataset': 'LUAD',
        'tip_figure': 'results/figures/tip_eu_vs_ricci.png',
        'sparsity_entropy_available': 'no',
        'support_size_available': 'no',
        'note': 'TIP is cached; per-bootstrap support vectors were not cached for post hoc sparsity/entropy recomputation.',
    },
    {
        'dataset': 'GSE161711',
        'tip_figure': 'results/figures/tip_eu_vs_ricci_cll_venetoclax_gse161711.png',
        'sparsity_entropy_available': 'no',
        'support_size_available': 'no',
        'note': 'Representative stability is available through TIP and H1 diagnostics; entropy/support-size recomputation would require rerunning topology loops.',
    },
    {
        'dataset': 'GSE165087 (bounded)',
        'tip_figure': 'results/figures/tip_eu_vs_ricci_cll_rs_scrna_gse165087.png',
        'sparsity_entropy_available': 'partial-only',
        'support_size_available': 'partial-only',
        'note': 'Bounded run preserves checkpoint and figure arrays, but no standalone entropy/support table was cached.',
    },
]
with (TABLES / 'representative_diagnostics_summary.csv').open('w', newline='') as f:
    writer = csv.DictWriter(
        f,
        fieldnames=['dataset', 'tip_figure', 'sparsity_entropy_available', 'support_size_available', 'note'],
    )
    writer.writeheader()
    writer.writerows(rep_rows)

labels = [row['dataset'] for row in summary_rows]
standard = np.array([row['standard_mean_f1'] for row in summary_rows], dtype=float)
harmonic = np.array([row['harmonic_mean_f1'] for row in summary_rows], dtype=float)
ccspr = np.array([row['ccspr_mean_f1'] for row in summary_rows], dtype=float)

x = np.arange(len(labels))
width = 0.24
fig, ax = plt.subplots(figsize=(8.5, 4.8))
ax.bar(x - width, standard, width, label='Standard', color='#9d6b53')
ax.bar(x, harmonic, width, label='Euclidean harmonic', color='#8da0cb')
ax.bar(x + width, ccspr, width, label='Ricci sparse', color='#66a61e')
ax.set_ylabel('Mean weighted F1')
ax.set_title('Implemented geometry study summary across completed datasets')
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=0)
ax.set_ylim(0, 1.05)
ax.legend(frameon=False, ncol=3, loc='upper center')
for xpos, val in zip(x, ccspr - harmonic):
    ax.text(xpos + width, ccspr[list(x).index(xpos)] + 0.025, f"Δ={val:+.3f}", ha='center', va='bottom', fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES / 'geometry_backend_summary.png', dpi=200, bbox_inches='tight')
plt.close(fig)

print('Wrote:')
for path in [
    TABLES / 'geometry_backend_summary.csv',
    TABLES / 'geometry_runtime_cost_summary.csv',
    TABLES / 'geometry_diagnostics_summary.csv',
    TABLES / 'representative_diagnostics_summary.csv',
    FIGURES / 'geometry_backend_summary.png',
]:
    print(path)

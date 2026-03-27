from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path('/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr')
RESULTS = ROOT / 'results'
TABLES = RESULTS / 'tables'
FIGURES = RESULTS / 'figures'
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

DATASETS = {
    'LUAD': RESULTS / 'ablation_summary.csv',
    'GSE161711': RESULTS / 'ablation_summary_cll_venetoclax_gse161711.csv',
}

frames = []
for dataset, path in DATASETS.items():
    df = pd.read_csv(path)
    df['dataset'] = dataset
    frames.append(df)
all_df = pd.concat(frames, ignore_index=True)
all_df.to_csv(TABLES / 'ablation_profiles_combined.csv', index=False)

best_rows = []
for dataset, df_ds in all_df.groupby('dataset'):
    for ablation, df_ab in df_ds.groupby('ablation'):
        row = df_ab.sort_values(['mean', 'n'], ascending=[False, False]).iloc[0]
        best_rows.append(
            {
                'dataset': dataset,
                'ablation': ablation,
                'best_value': row['ablation_value'],
                'best_mean_f1': row['mean'],
                'std': row['std'],
                'n': int(row['n']),
            }
        )
best_df = pd.DataFrame(best_rows)
best_df.to_csv(TABLES / 'ablation_best_summary.csv', index=False)

# Separate tables per dataset for easy manuscript inclusion.
for dataset in DATASETS:
    best_df[best_df['dataset'] == dataset].to_csv(
        TABLES / f"ablation_best_{dataset.lower().replace('(', '').replace(')', '').replace(' ', '_')}.csv",
        index=False,
    )

ablations = ['distance_mode', 'iters', 'k', 'lambda', 'top_k', 'normalize']
fig, axes = plt.subplots(3, 2, figsize=(11, 11))
axes = axes.flatten()
colors = {'LUAD': '#8da0cb', 'GSE161711': '#66a61e'}
for ax, ablation in zip(axes, ablations):
    subset = all_df[all_df['ablation'] == ablation].copy()
    for dataset, df_ds in subset.groupby('dataset'):
        x = df_ds['ablation_value'].astype(str).tolist()
        y = df_ds['mean'].astype(float).tolist()
        ax.plot(range(len(x)), y, marker='o', linewidth=2, label=dataset, color=colors.get(dataset))
        ax.set_xticks(range(len(x)))
        ax.set_xticklabels(x, rotation=30, ha='right')
    ax.set_title(ablation)
    ax.set_ylabel('Mean weighted F1')
    ax.grid(alpha=0.25)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
fig.suptitle('Completed LUAD and GSE161711 ablation profiles (reused, not rerun)', y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGURES / 'ablation_profiles_luad_gse161711.png', dpi=200, bbox_inches='tight')
plt.close(fig)

# One compact bar chart for best-value gains by ablation family.
fig, ax = plt.subplots(figsize=(10, 4.5))
plot_df = best_df.copy()
plot_df['label'] = plot_df['dataset'] + ':' + plot_df['ablation']
ax.bar(plot_df['label'], plot_df['best_mean_f1'], color=['#8da0cb' if d == 'LUAD' else '#66a61e' for d in plot_df['dataset']])
ax.set_ylabel('Best mean weighted F1')
ax.set_title('Best completed ablation values by dataset and family')
ax.set_xticklabels(plot_df['label'], rotation=45, ha='right')
fig.tight_layout()
fig.savefig(FIGURES / 'ablation_best_summary.png', dpi=200, bbox_inches='tight')
plt.close(fig)

print('Wrote ablation artifact tables and figures')

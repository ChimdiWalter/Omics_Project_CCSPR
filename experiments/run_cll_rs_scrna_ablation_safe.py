#!/usr/bin/env python3
from __future__ import annotations

import argparse
import resource
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score
from sklearn.model_selection import StratifiedShuffleSplit

from ccspr.datasets.cll_rs_scrna import load_cll_rs_scrna
from ccspr.eval.classification import FeatureParams, _evaluate_feature_model, _fit_predict_logreg
from ccspr.plots.figures import plot_f1_bars_with_ci, plot_lambda_robustness, plot_lifetime_prominence, plot_tip_eu_vs_ricci
from ccspr.preprocess.basic import select_top_variable_genes, standardize
from ccspr.utils.io import ensure_dir, read_yaml, save_json


def _ts() -> str:
    return datetime.now().isoformat(timespec='seconds')


def _rss_gb() -> float:
    return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / (1024.0 * 1024.0)


def _log(msg: str) -> None:
    print(f'[{_ts()}] {msg}', flush=True)


def _summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=['model', 'mean', 'std', 'n', 'ci95'])
    out = (
        df.groupby('model')['f1_weighted']
        .agg(['mean', 'std', 'count'])
        .reset_index()
        .rename(columns={'count': 'n'})
    )
    out['ci95'] = 1.96 * out['std'].fillna(0.0) / np.sqrt(out['n'].clip(lower=1))
    return out


def _ablation_summary(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=['ablation', 'ablation_value', 'mean', 'std', 'n'])
    return (
        df.groupby(['ablation', 'ablation_value'])['f1_weighted']
        .agg(['mean', 'std', 'count'])
        .reset_index()
        .rename(columns={'count': 'n'})
    )


def _write_checkpoint(path: Path, phase: str, meta: dict) -> None:
    payload = {'timestamp': _ts(), 'phase': phase, 'max_rss_gb': round(_rss_gb(), 3)}
    payload.update(meta)
    save_json(payload, path)


def _write_partial(out_dir: Path, tables_dir: Path, tag: str, metrics_main: pd.DataFrame, ablation_df: pd.DataFrame, aux_state: dict[str, np.ndarray], figure_npz: Path) -> None:
    metrics_main.to_csv(out_dir / f'metrics_{tag}.partial.csv', index=False)
    _summary(metrics_main).to_csv(out_dir / f'main_results_{tag}.partial.csv', index=False)
    _summary(metrics_main).to_csv(tables_dir / f'main_results_{tag}.partial.csv', index=False)
    ab_sum = _ablation_summary(ablation_df)
    ab_sum.to_csv(out_dir / f'ablation_summary_{tag}.partial.csv', index=False)
    ab_sum.to_csv(tables_dir / f'ablation_summary_{tag}.partial.csv', index=False)
    np.savez(
        figure_npz,
        tip_eu=np.asarray(aux_state.get('tip_eu', np.array([], dtype=float)), dtype=float),
        tip_ricci=np.asarray(aux_state.get('tip_ricci', np.array([], dtype=float)), dtype=float),
        lifetime_eu=np.asarray(aux_state.get('lifetime_eu', np.array([], dtype=float)), dtype=float),
        lifetime_ricci=np.asarray(aux_state.get('lifetime_ricci', np.array([], dtype=float)), dtype=float),
        prominence_eu=np.asarray(aux_state.get('prominence_eu', np.array([], dtype=float)), dtype=float),
        prominence_ricci=np.asarray(aux_state.get('prominence_ricci', np.array([], dtype=float)), dtype=float),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='configs/cll_rs_scrna_ablation_safe.yaml')
    args = ap.parse_args()
    cfg = read_yaml(args.config)

    seed = int(cfg.get('seed', 42))
    np.random.seed(seed)

    out = Path(cfg['output'].get('results_dir', 'results'))
    tag = str(cfg['output'].get('tag', 'cll_rs_scrna_gse165087_ablation_safe'))
    tables = out / 'tables'
    figures = out / 'figures'
    ensure_dir(out)
    ensure_dir(tables)
    ensure_dir(figures)
    ensure_dir(Path(cfg.get('cache_dir', 'data/cache')))

    final_metrics = out / f'metrics_{tag}.csv'
    final_main = out / f'main_results_{tag}.csv'
    final_ablation = out / f'ablation_summary_{tag}.csv'
    if bool(cfg.get('skip_if_complete', True)) and final_metrics.exists() and final_main.exists() and final_ablation.exists():
        _log(f'Final outputs already exist for {tag}; skipping ablation-safe run')
        return

    checkpoint_json = out / f'checkpoint_{tag}.json'
    dataset_meta_json = out / f'dataset_metadata_{tag}.json'
    figure_npz = out / cfg.get('figure_arrays', {}).get('filename', f'figure_arrays_{tag}.npz')
    save_json({'timestamp': datetime.utcnow().isoformat() + 'Z', 'config': cfg}, out / f'run_{tag}.json')
    _write_checkpoint(checkpoint_json, 'starting', {'config': args.config})

    _log('Loading bounded scRNA dataset for ablation-safe pass')
    ds = load_cll_rs_scrna(
        data_path=cfg['dataset']['data_path'],
        label_key=cfg['dataset'].get('label_key', 'label'),
        min_genes=int(cfg['dataset'].get('min_genes', 200)),
        min_cells=int(cfg['dataset'].get('min_cells', 3)),
        n_hvg=int(cfg['dataset'].get('n_hvg', 1000)),
        n_pcs=int(cfg['dataset'].get('n_pcs', 25)),
        max_cells=cfg['dataset'].get('max_cells', 400),
        seed=seed,
        use_preprocessed_if_available=bool(cfg['dataset'].get('use_preprocessed_if_available', True)),
    )
    save_json(ds['meta'], dataset_meta_json)
    _write_checkpoint(checkpoint_json, 'dataset_loaded', {'n_cells': int(ds['X'].shape[0]), 'n_features': int(ds['X'].shape[1])})
    _log(f"Dataset ready: cells={ds['X'].shape[0]} features={ds['X'].shape[1]} max_rss_gb={_rss_gb():.2f}")

    harmonic_params = FeatureParams(**cfg['protocol']['harmonic_params'])
    ccspr_params = FeatureParams(**cfg['protocol']['ccspr_params'])
    var_genes = min(int(cfg['protocol'].get('var_genes', ds['X'].shape[1])), ds['X'].shape[1])
    splitter = StratifiedShuffleSplit(n_splits=int(cfg['protocol'].get('n_splits', 1)), test_size=float(cfg['protocol'].get('test_size', 0.25)), random_state=seed)
    itr, ite = next(splitter.split(ds['X'], ds['y']))
    x_train = ds['X'][itr]
    x_test = ds['X'][ite]
    y_train = ds['y'][itr]
    y_test = ds['y'][ite]
    x_train_var, idx_var = select_top_variable_genes(x_train, top_n=var_genes)
    x_test_var = x_test[:, idx_var]

    metrics_main = pd.DataFrame()
    ablation_rows = []
    aux_state: dict[str, np.ndarray] = {}

    _log('Running standard baseline')
    xtr_std, xts_std = standardize(x_train_var, x_test_var)
    y_pred_std = _fit_predict_logreg(xtr_std, y_train, xts_std, seed)
    metrics_main = pd.concat([metrics_main, pd.DataFrame([{'split': 0, 'model': 'standard', 'f1_weighted': float(f1_score(y_test, y_pred_std, average='weighted')), 'distance_mode': 'euclidean', 'solver': 'none'}])], ignore_index=True)
    _write_partial(out, tables, tag, metrics_main, pd.DataFrame(ablation_rows), aux_state, figure_npz)

    _log('Running harmonic baseline')
    f1_h, out_h = _evaluate_feature_model(x_train_var, y_train, x_test_var, y_test, harmonic_params, seed + 1000, cfg.get('cache_dir', 'data/cache'))
    aux_state['tip_eu'] = out_h['tip']
    aux_state['lifetime_eu'] = out_h['lifetime']
    aux_state['prominence_eu'] = out_h['prominence']
    metrics_main = pd.concat([metrics_main, pd.DataFrame([{'split': 0, 'model': 'harmonic', 'f1_weighted': float(f1_h), 'distance_mode': harmonic_params.mode, 'solver': harmonic_params.solver, 'k': harmonic_params.k, 'iters': harmonic_params.iters, 'lambda': harmonic_params.lambda_, 'top_k': harmonic_params.top_k, 'normalize': harmonic_params.normalize}])], ignore_index=True)
    _write_partial(out, tables, tag, metrics_main, pd.DataFrame(ablation_rows), aux_state, figure_npz)

    _log('Running CC-SPR main comparison')
    f1_c, out_c = _evaluate_feature_model(x_train_var, y_train, x_test_var, y_test, ccspr_params, seed + 2000, cfg.get('cache_dir', 'data/cache'))
    aux_state['tip_ricci'] = out_c['tip']
    aux_state['lifetime_ricci'] = out_c['lifetime']
    aux_state['prominence_ricci'] = out_c['prominence']
    metrics_main = pd.concat([metrics_main, pd.DataFrame([{'split': 0, 'model': 'ccspr', 'f1_weighted': float(f1_c), 'distance_mode': ccspr_params.mode, 'solver': ccspr_params.solver, 'k': ccspr_params.k, 'iters': ccspr_params.iters, 'lambda': ccspr_params.lambda_, 'top_k': ccspr_params.top_k, 'normalize': ccspr_params.normalize}])], ignore_index=True)
    _write_partial(out, tables, tag, metrics_main, pd.DataFrame(ablation_rows), aux_state, figure_npz)

    if bool(cfg.get('ablation', {}).get('enabled', False)):
        kind = str(cfg['ablation'].get('kind', 'lambda'))
        values = list(cfg['ablation'].get('values', []))
        for idx, value in enumerate(values):
            params = FeatureParams(**{**cfg['protocol']['ccspr_params'], 'lambda_': float(value)}) if kind == 'lambda' else FeatureParams(**{**cfg['protocol']['ccspr_params'], 'k': int(value)})
            _log(f'Running bounded scRNA ablation {kind}={value}')
            f1_a, _ = _evaluate_feature_model(x_train_var, y_train, x_test_var, y_test, params, seed + 3000 + idx, cfg.get('cache_dir', 'data/cache'))
            ablation_rows.append({'split': 0, 'model': 'ccspr', 'f1_weighted': float(f1_a), 'distance_mode': params.mode, 'solver': params.solver, 'k': params.k, 'iters': params.iters, 'lambda': params.lambda_, 'top_k': params.top_k, 'normalize': params.normalize, 'ablation': kind, 'ablation_value': value, 'phase': 'ablation', 'dataset': tag})
            _write_partial(out, tables, tag, metrics_main, pd.DataFrame(ablation_rows), aux_state, figure_npz)
            _write_checkpoint(checkpoint_json, f'ablation_{kind}_{value}', {'rows_written': len(ablation_rows)})
            _log(f'Completed bounded scRNA ablation {kind}={value}; max_rss_gb={_rss_gb():.2f}')

    metrics_main = metrics_main.copy()
    metrics_main['phase'] = 'main'
    metrics_main['dataset'] = tag
    metrics_main.to_csv(final_metrics, index=False)
    summary = _summary(metrics_main)
    summary.to_csv(final_main, index=False)
    summary.to_csv(tables / f'main_results_{tag}.csv', index=False)

    ablation_df = pd.DataFrame(ablation_rows)
    ablation_summary = _ablation_summary(ablation_df)
    ablation_summary.to_csv(final_ablation, index=False)
    ablation_summary.to_csv(tables / f'ablation_summary_{tag}.csv', index=False)

    if aux_state.get('tip_eu', np.array([])).size and aux_state.get('tip_ricci', np.array([])).size:
        plot_tip_eu_vs_ricci(aux_state['tip_eu'], aux_state['tip_ricci'], str(figures / f'tip_eu_vs_ricci_{tag}.png'))
    if aux_state.get('lifetime_eu', np.array([])).size and aux_state.get('lifetime_ricci', np.array([])).size:
        plot_lifetime_prominence(aux_state['lifetime_eu'], aux_state['lifetime_ricci'], aux_state.get('prominence_eu', np.array([], dtype=float)), aux_state.get('prominence_ricci', np.array([], dtype=float)), str(figures / f'h1_lifetime_eu_vs_ricci_{tag}.png'), str(figures / f'h1_prominence_eu_vs_ricci_{tag}.png'))
    plot_f1_bars_with_ci(metrics_main, str(figures / f'f1_bars_ci_{tag}.png'))
    if not ablation_df.empty and str(cfg['ablation'].get('kind', 'lambda')) == 'lambda':
        plot_lambda_robustness(ablation_df, str(figures / f'lambda_robustness_{tag}.png'))

    _write_checkpoint(checkpoint_json, 'complete', {'final_metrics': str(final_metrics), 'final_main_results': str(final_main), 'final_ablation_summary': str(final_ablation)})
    _log(f'Completed scRNA ablation-safe run; outputs at {final_metrics}')


if __name__ == '__main__':
    main()

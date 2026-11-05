"""Rebuild aggregate HOI-TRUST figures without patient-level data or model fitting.

Figure 3 uses a documented coarse grid with a minimum retained-count increment.
These are aggregate-only reproductions; fine-grained private curve points are omitted.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, MaxNLocator

COLORS = {'patient': '#df7900', 'environmental': '#2462a5'}
ORDER = ['patient', 'environmental']
PATHOGENS = {
    'DS_E_coli': 'Drug-susceptible E. coli', 'ESBL_E_coli': 'ESBL E. coli',
    'DS_K_pneumoniae': 'Drug-susceptible K. pneumoniae',
    'ESBL_K_pneumoniae': 'ESBL K. pneumoniae', 'C_diff': 'C. difficile',
    'VSE_faecalis': 'Vancomycin-susceptible E. faecalis',
    'VRE_faecium': 'Vancomycin-resistant E. faecium', 'MSSA': 'MSSA', 'MRSA': 'MRSA',
    'DS_P_aeruginosa': 'Drug-susceptible P. aeruginosa',
    'DR_P_aeruginosa': 'Drug-resistant P. aeruginosa',
}

def render_all(root):
    root = Path(root)
    tables, figures = root / 'results', root / 'figures'
    figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.dpi': 300})

    def read(name):
        return pd.read_csv(tables / name)

    def save(fig, name):
        fig.savefig(figures / (name + '.png'), bbox_inches='tight', facecolor='white')
        fig.savefig(figures / (name + '.pdf'), bbox_inches='tight', facecolor='white')
        plt.close(fig)

    # Primary component-bootstrap estimates and paired differences.
    absolute = read('component_bootstrap_95ci.csv')
    paired = read('paired_stack_vs_best_single_95ci.csv')
    fig, ax = plt.subplots(1, 3, figsize=(16, 6), layout='constrained',
                           gridspec_kw={'width_ratios': [1, 1, 1.3]})
    combinations = [('patient', 'Strict nested stack'), ('patient', 'Best single model: CatBoost'),
                    ('environmental', 'Strict nested stack'), ('environmental', 'Best single model: Extra Trees')]
    labels = ['Patient — stack', 'Patient — CatBoost', 'Environmental — stack', 'Environmental — Extra Trees']
    for j, metric in enumerate(['AUROC', 'AUPRC']):
        for y, (pathway, model) in enumerate(combinations):
            row = absolute.loc[(absolute.pathway == pathway) & (absolute.model == model)
                               & (absolute.metric == metric)].iloc[0]
            x, lo, hi = row['estimate'], row.ci_95_lower, row.ci_95_upper
            ax[j].errorbar(x, y, xerr=[[x-lo], [hi-x]], fmt='o' if 'stack' in model else 's',
                           color=COLORS[pathway], capsize=3)
            ax[j].annotate(f'{x:.4f}', (x, y), xytext=(6, -12), textcoords='offset points', fontsize=9)
        ax[j].set(title=metric, xlim=(0, 1), yticks=range(4))
        ax[j].set_yticklabels(labels if j == 0 else [])
        ax[j].invert_yaxis()
        ax[j].set_xlabel('Estimate (95% component-bootstrap interval)')
        ax[j].grid(axis='x', alpha=.15)
    metrics = ['AUROC', 'AUPRC', 'precision', 'recall', 'enrichment']
    for pathway, shift in [('patient', -.11), ('environmental', .11)]:
        for y, metric in enumerate(metrics):
            row = paired.loc[(paired.pathway == pathway) & (paired.metric == metric)].iloc[0]
            x, lo, hi = row.paired_difference, row.paired_ci_95_lower, row.paired_ci_95_upper
            ax[2].errorbar(x, y+shift, xerr=[[x-lo], [hi-x]], fmt='o', color=COLORS[pathway],
                           capsize=3, label=pathway.title() if y == 0 else None)
    ax[2].axvline(0, color='.4', ls='--')
    ax[2].set(title='Paired differences', yticks=range(5),
              yticklabels=['AUROC', 'AUPRC', 'Precision at 10%', 'Recall at 10%', 'Enrichment at 10%'],
              xlabel='Stack minus pathway-specific comparator')
    ax[2].invert_yaxis(); ax[2].legend(frameon=False)
    fig.suptitle('Primary internal validation: estimates and paired differences')
    save(fig, 'results_primary_component_bootstrap')

    # Single-repeat sensitivities; uncertainty intervals were not estimated.
    data = read('strict_sensitivity_results.csv')
    data = data.loc[data.model == 'Strict nested stack']
    scenarios = ['one_exact_copy', 'adults_only', 'exclude_irregular_matched_sets',
                 'exclude_resolved_any_abx_summary', 'exclude_largest_dependency_component']
    names = ['One exact copy', 'Adults only', 'Exclude irregular matched sets',
             'Remove antibiotic summary (patient only)', 'Exclude largest component']
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), sharey=True, layout='constrained')
    for axis, key, title in zip(axes, ['delta_vs_primary_AUROC', 'delta_vs_primary_AUPRC',
                                     'delta_vs_primary_precision'], ['Change in AUROC', 'Change in AUPRC',
                                                                     'Change in precision at 10%']):
        for pathway, shift in [('patient', -.10), ('environmental', .10)]:
            sub = data.loc[data.pathway == pathway].set_index('scenario')
            for y, scenario in enumerate(scenarios):
                if scenario in sub.index:
                    axis.scatter(sub.loc[scenario, key], y+shift, color=COLORS[pathway],
                                 label=pathway.title() if y == 0 else None)
        axis.axvline(0, color='.4', ls='--'); axis.set_title(title)
        axis.set_xlabel('Difference from corresponding primary repeat')
        axis.xaxis.set_major_locator(MaxNLocator(4))
        axis.grid(axis='x', alpha=.15)
    axes[0].set(yticks=range(5), yticklabels=names); axes[0].invert_yaxis()
    axes[-1].legend(frameon=False)
    fig.suptitle('Development sensitivities: single-repeat diagnostics without confidence intervals')
    save(fig, 'results_sensitivity_analyses')

    # Capacity points and coarsened selective error curve, with no per-record thresholds.
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for pathway in ORDER:
        cap = read(f'capacity_curve_operating_{pathway}.csv').sort_values('capacity_fraction')
        curve = read(f'selective_curve_public_{pathway}.csv').sort_values('coverage')
        axes[0].plot(cap.capacity_fraction, cap.recall_at_capacity, '-o', color=COLORS[pathway], label=pathway.title())
        axes[1].plot(curve.coverage, curve.selective_error, color=COLORS[pathway], label=pathway.title())
    axes[0].axhline(.10, color='.4', ls='--'); axes[0].axvline(.10, color='.4', ls=':')
    axes[1].axhline(.20, color='.4', ls='--')
    for axis in axes:
        axis.xaxis.set_major_formatter(PercentFormatter(1))
        axis.yaxis.set_major_formatter(PercentFormatter(1)); axis.legend(frameon=False)
        axis.grid(alpha=.15)
    axes[0].set(title='Case capture by review capacity', xlabel='Nominal within-pathogen review fraction', ylabel='Case recall')
    axes[1].set(title='Selective error on the public aggregate grid', xlabel='Retained fraction', ylabel='Classification error')
    fig.suptitle('Operating-policy development; the 20% capacity point is descriptive')
    save(fig, 'operating_policy_curves')

    # Component and observation coverage and set composition.
    uncertainty = read('v2_legacy_uncertainty_metrics_NOT_CONFIRMATORY.csv').set_index('pathway').loc[ORDER]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
    x = np.arange(2)
    for shift, key, label, color in [(-.16, 'observation_coverage', 'Observation', '#4c78a8'),
                                     (.16, 'component_complete_coverage', 'Complete component', '#df7900')]:
        bars = axes[0].bar(x+shift, uncertainty[key], .30, label=label, color=color)
        axes[0].bar_label(bars, labels=[f'{v:.4f}' for v in uncertainty[key]], padding=3)
    axes[0].axhline(.90, ls='--', color='.4', label='Nominal 90%')
    bottom = np.zeros(2)
    for key, label, color in [('singleton_fraction', 'Singleton', '#4c78a8'),
                               ('ambiguous_fraction', 'Ambiguous', '#f2a43c'), ('empty_fraction', 'Empty', '#e45756')]:
        v = uncertainty[key].to_numpy()
        axes[1].bar(x, v, bottom=bottom, label=label, color=color)
        for j in range(2):
            if v[j] > .001: axes[1].text(j, bottom[j]+v[j]/2, f'{v[j]:.2%}', ha='center', va='center')
        bottom += v
    for axis, title in zip(axes, ['Empirical coverage', 'Prediction-set composition']):
        axis.set(xticks=x, xticklabels=[p.title() for p in ORDER], title=title)
        axis.yaxis.set_major_formatter(PercentFormatter(1))
        axis.legend(frameon=False, loc='upper center', bbox_to_anchor=(.5, -.10), ncol=3, fontsize=9)
    axes[0].set_ylim(0, 1.09); axes[1].set_ylim(0, 1)
    fig.suptitle('Conformal uncertainty in the non-confirmatory compatibility partition')
    save(fig, 'results_conformal_compatibility')

    # Pathogen-specific lift.
    pathogens = read('legacy_within_pathogen_metrics_NOT_CONFIRMATORY.csv')
    fig, axis = plt.subplots(figsize=(12, 8), layout='constrained')
    for pathway, shift in [('patient', -.10), ('environmental', .10)]:
        sub = pathogens.loc[pathogens.pathway == pathway].set_index('run')
        for y, code in enumerate(PATHOGENS):
            row = sub.loc[code]
            axis.scatter(row.AUPRC_lift, y+shift, s=max(20, float(row.observations)/12),
                         color=COLORS[pathway], alpha=.85, label=pathway.title() if y == 0 else None)
    axis.axvline(1, color='.4', ls='--')
    axis.set(yticks=range(11), yticklabels=list(PATHOGENS.values()),
             xlabel='AUPRC / pathogen-specific sampled case proportion',
             title='Within-pathogen AUPRC lift: non-confirmatory compatibility partition')
    axis.invert_yaxis(); axis.legend(frameon=False)
    save(fig, 'results_within_pathogen_heterogeneity')

    # Grouped permutation percentiles are not confidence intervals.
    reliance = read('v2_locked_candidate_grouped_permutation_importance.csv')
    domains = ['Colonisation pressure', 'Pathogen context and interactions', 'Comorbidity',
               'Demographic and surgery', 'Prior antibiotic exposure']
    fig, axis = plt.subplots(figsize=(12, 6), layout='constrained')
    for pathway, shift in [('patient', -.12), ('environmental', .12)]:
        sub = reliance.loc[reliance.pathway == pathway].set_index('domain')
        labelled = False
        for y, domain in enumerate(domains):
            if domain not in sub.index: continue
            row = sub.loc[domain]; x = row.mean_AUPRC_decrease
            axis.errorbar(x, y+shift, xerr=[[x-row.lower_2_5], [row.upper_97_5-x]], fmt='o',
                          color=COLORS[pathway], capsize=4, label=pathway.title() if not labelled else None)
            labelled = True
    axis.axvline(0, color='.4', ls='--')
    axis.set(yticks=range(5), yticklabels=domains, xlabel='AUPRC decrease after joint domain permutation',
             title='Grouped permutation reliance on the operating partition')
    axis.invert_yaxis(); axis.legend(frameon=False)
    fig.supxlabel('Bars are empirical permutation percentiles, not confidence intervals. Associations are not causal.', fontsize=9)
    save(fig, 'results_grouped_permutation_reliance')

    # Dashboard figure uses two aggregate margins, without individual scores.
    queue = read('dashboard_review_queue.csv')
    composition = read('dashboard_uncertainty_composition.csv')
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), layout='constrained')
    for pathway, shift in [('patient', -.16), ('environmental', .16)]:
        sub = queue.loc[queue.pathway == pathway].set_index('pathogen_code')
        vals = [sub.loc[code, 'selected_for_review'] for code in PATHOGENS]
        axes[0].barh(np.arange(11)+shift, vals, .30, color=COLORS[pathway], label=pathway.title())
    axes[0].set(yticks=range(11), yticklabels=list(PATHOGENS.values()), xlabel='Selected observations', title='Review queue by pathogen target')
    axes[0].invert_yaxis(); axes[0].legend(frameon=False)
    bottom = np.zeros(2)
    totals = composition.groupby('pathway').observations.sum()
    for label, color in [('Singleton', '#4c78a8'), ('Ambiguous', '#f2a43c'), ('Empty', '#e45756')]:
        sub = composition.loc[composition.uncertainty == label].set_index('pathway')
        vals = np.array([sub.loc[p, 'observations']/totals[p] if p in sub.index else 0 for p in ORDER])
        axes[1].bar(np.arange(2), vals, bottom=bottom, label=label, color=color); bottom += vals
    axes[1].set(xticks=range(2), xticklabels=[p.title() for p in ORDER], title='Uncertainty composition', ylim=(0, 1))
    axes[1].yaxis.set_major_formatter(PercentFormatter(1)); axes[1].legend(frameon=False)
    selected = int(queue.selected_for_review.sum()); n = int(queue.observations.sum())
    ambiguous = int(composition.loc[composition.uncertainty == 'Ambiguous', 'observations'].sum())
    fig.suptitle(f'Retrospective surveillance review — {n:,} observations; {selected:,} selected ({selected/n:.1%}); {ambiguous/n:.1%} ambiguous')
    save(fig, 'results_dashboard_demonstrator')
    print(f'Saved seven aggregate figure reproductions to {figures}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    render_all(args.root)

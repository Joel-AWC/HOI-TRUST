"""Package the completed HOI-TRUST run for a public GitHub/Zenodo release.

Run this file in a SEPARATE Colab notebook. It reads saved files, never fits models,
and never uploads anything. A missing file or a failed run check stops the export.
The output ZIP is a release candidate for inspection before publication.
"""
from pathlib import Path
from datetime import datetime, timezone
import ast
import csv
import hashlib
import io
import json
import math
import os
import shutil
import tempfile
import zipfile

# Change BASE_DIR only if your saved project is in a different Drive folder.
BASE_DIR = Path(os.environ.get('HOI_TRUST_DATA_DIR',
    '/content/drive/MyDrive/MSc_Project_Private_Data'))
# Optional: set an exact notebook path if its filename differs from these defaults.
NOTEBOOK_FILE = None
OUTPUT_PARENT = None  # Defaults to a sibling folder named HOI_TRUST_public_release.
RELEASE_VERSION = 'v1.0.1'
DOWNLOAD_ZIP = True

EXPORTER_MARKER = 'HOI_TRUST_AGGREGATE_RELEASE_EXPORTER_1'
EXPECTED_FINGERPRINT = '58ad1bed6db181516c44dd2a405e24017d983d2c3090906780ee00c9a9f4725a'
EXPECTED_CODE_VERSION = '89a2d6e1875e10913bbdcd03f04c32c155034a4618b992e51fb236f51c27c8df'
EXPECTED_CELL_SOURCE_SHA = '7df61d8c649998d84cff3bce6c1218fe290431dff736c81ce9c77b5c699ea927'
EXPECTED_DIRECT_SHA = '0e20a5a7ffb858580ec9cd09a5d684be7f9376329b10330c0c8a304548df75f5'
EXPECTED_ENV_SHA = '4f13bd59acf7cffaa2194c6f124e0ed6a15f627efb044e9dc9d9ac09aa0056d2'
EXPECTED_PYTHON = '3.13.16'
MIN_CURVE_BLOCK = 30
PATHWAYS = {'patient', 'environmental'}
CANDIDATES = {'Strict nested stack', 'CatBoost', 'Extra Trees', 'Explainable Boosting Machine',
              'XGBoost', 'Logistic regression', 'Pairwise ranker', 'Multi-task neural network'}
PATHOGENS = {'DS_E_coli', 'ESBL_E_coli', 'DS_K_pneumoniae', 'ESBL_K_pneumoniae',
             'C_diff', 'VSE_faecalis', 'VRE_faecium', 'MSSA', 'MRSA', 'DS_P_aeruginosa', 'DR_P_aeruginosa'}
FORBIDDEN = {'patient_id', 'deidentified_patient_id', 'pathway_patient_id', 'matched_set_id',
             'dependency_component', 'group_index', 'group_binary', 'record_id', 'y_true',
             'y_pred', 'raw_score', 'ranking_score', 'model_score', 'v2_raw_score',
             'v2_calibrated_score', 'subject_id', 'hadm_id'}
METRICS = ('AUROC AUPRC AUPRC_lift Brier Brier_skill ECE_10 calibration_intercept '
           'calibration_slope prevalence').split()
CAPACITY = ('capacity_fraction nominal_n actual_n_including_ties true_positives '
            'precision_at_capacity recall_at_capacity enrichment_at_capacity '
            'number_needed_to_review work_saved_vs_review_all selected_fraction').split()

class ExportError(RuntimeError):
    pass

def require(ok, message):
    if not ok:
        raise ExportError(message)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def number(row, key):
    try:
        value = float(row[key])
    except (KeyError, ValueError, TypeError):
        raise ExportError(f'Missing or non-numeric aggregate field: {key}')
    require(math.isfinite(value), f'Non-finite required aggregate field: {key}')
    return value

def source_bytes(path, limit=40_000_000):
    require(path.is_file() and not path.is_symlink(), f'Missing ordinary saved file: {path.name}')
    require(path.stat().st_size <= limit, f'Unexpected file size: {path.name}')
    return path.read_bytes()

def rows_from(path, required, maximum=250, keep=None):
    raw = source_bytes(path)
    text = raw.decode('utf-8-sig')
    reader = csv.DictReader(io.StringIO(text))
    columns = reader.fieldnames or []
    require(len(columns) == len(set(columns)), f'Duplicate CSV headings: {path.name}')
    require(not (set(c.lower() for c in columns) & FORBIDDEN),
            f'Record-level column detected; file excluded: {path.name}')
    require(set(required) <= set(columns), f'Unexpected schema in {path.name}; required: {", ".join(required)}')
    rows = list(reader)
    require(0 < len(rows) <= maximum, f'Unexpected aggregate row count in {path.name}: {len(rows)}')
    require(all(None not in row for row in rows), f'Malformed CSV rows: {path.name}')
    for row in rows:
        if 'pathway' in row:
            row['pathway'] = row['pathway'].lower()
            require(row['pathway'] in PATHWAYS, f'Unexpected pathway in {path.name}')
        if 'analysis' in row:
            require(row['analysis'] in PATHWAYS, f'Unexpected analysis branch in {path.name}')
        if 'candidate' in row:
            require(row['candidate'] in CANDIDATES, f'Unexpected candidate in {path.name}')
        if 'run' in row:
            require(row['run'] in PATHOGENS, f'Unexpected pathogen code in {path.name}')
    keep = list(dict.fromkeys(keep or required))
    keep = [key for key in keep if key in columns]
    require(set(required) <= set(keep), f'Exporter schema omits required fields: {path.name}')
    public = [{key: row[key] for key in keep} for row in rows]
    return public, {'source_file': path.name, 'source_sha256': sha(raw), 'source_rows': len(rows),
                    'exported_columns': keep, 'omitted_columns': [c for c in columns if c not in keep]}

def csv_bytes(rows, columns=None):
    require(bool(rows), 'An aggregate export would be empty.')
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=columns or list(rows[0]), lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return stream.getvalue().encode('utf-8')

def one(rows, **keys):
    found = [row for row in rows if all(row.get(k) == v for k, v in keys.items())]
    require(len(found) == 1, f'Expected one aggregate row for {keys}; found {len(found)}')
    return found[0]

def rounded_check(rows, column, expected, **keys):
    value = number(one(rows, **keys), column)
    require(abs(value - expected) <= .00005001,
            f'Manuscript mismatch for {keys}, {column}: observed {value}; expected approximately {expected}')

def snapshot_notebook(base, exact_path=None):
    candidates = [Path(exact_path)] if exact_path else [
        base/'HOI_TRUST_submission_ready.ipynb', base/'HOI_TRUST_submission_ready (5).ipynb',
        base/'HOI_TRUST_pipeline.ipynb']
    for path in candidates:
        if not path.is_file():
            continue
        raw = source_bytes(path)
        notebook = json.loads(raw)
        sources = [''.join(c.get('source', [])) for c in notebook.get('cells', []) if c.get('cell_type') == 'code']
        if len(sources) != 52 or sha('\n'.join(sources).encode()) != EXPECTED_CELL_SOURCE_SHA:
            continue
        # Reconstruct the direct lock from literal pins, without executing notebook code.
        pins = None
        for source in sources:
            if 'PINNED_REQUIREMENTS = {' not in source:
                continue
            for node in ast.parse(source).body:
                if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'PINNED_REQUIREMENTS' for t in node.targets):
                    pins = ast.literal_eval(node.value)
        require(isinstance(pins, dict), 'Cannot read the recorded direct package pins.')
        lock = ('\n'.join(f'{k}=={v}' for k, v in sorted(pins.items()))+'\n').encode()
        require(sha(lock) == EXPECTED_DIRECT_SHA, 'The notebook package pins differ from the manuscript run.')
        # Keep source only. Strip all cell metadata, outputs, execution counts and attachments.
        clean_cells = []
        for cell in notebook['cells']:
            if cell['cell_type'] not in {'code', 'markdown'}:
                continue
            clean = {'cell_type': cell['cell_type'], 'metadata': {}, 'source': cell.get('source', [])}
            if clean['cell_type'] == 'code':
                clean.update(outputs=[], execution_count=None)
            clean_cells.append(clean)
        clean_cells.insert(1, {'cell_type': 'markdown', 'metadata': {}, 'source': [
            '## Public release\n',
            'Saved outputs, execution counts, attachments and widget state have been removed.\n',
            'Code-cell sources match the completed manuscript run. Obtain authorised PhysioNet access before execution.\n']})
        public = {'cells': clean_cells, 'nbformat': 4, 'nbformat_minor': 5,
                  'metadata': {k: notebook.get('metadata', {})[k] for k in ['kernelspec', 'language_info'] if k in notebook.get('metadata', {})}}
        return json.dumps(public, ensure_ascii=False, indent=1).encode()+b'\n', lock, {
            'source_notebook_file': path.name, 'source_notebook_sha256': sha(raw),
            'code_cells': len(sources), 'code_cell_source_sha256': EXPECTED_CELL_SOURCE_SHA}
    raise ExportError('No matching saved analysis notebook was found. Set NOTEBOOK_FILE to the retained notebook. '
                      'Do not append this exporter to that notebook or change the expected run hashes.')

def coarsen_curve(rows, policy_row=None):
    # Each public point comes from a saved aggregate curve; no interpolation or label reconstruction.
    by_n = {}
    totals = []
    for row in rows:
        n = number(row, 'selected_n'); coverage = number(row, 'coverage'); error = number(row, 'selective_error')
        require(n.is_integer() and 0 < coverage <= 1 and 0 <= error <= 1, 'Invalid saved selective-curve values.')
        n = int(n); totals.append(n/coverage)
        item = {'selected_n': n, 'coverage': coverage, 'selective_error': error}
        if n in by_n:
            require(abs(by_n[n]['selective_error']-error) < 1e-10, 'Inconsistent selective errors at the same retained count.')
        by_n[n] = item
    total = round(max(totals))
    require(all(abs(v-total) < .01 for v in totals), 'Selective curve does not use a single operating-partition denominator.')
    require(total in {4874, 4346}, 'Unexpected selective operating-partition size.')
    possible = [r for r in by_n.values() if r['selected_n'] >= MIN_CURVE_BLOCK]
    possible.sort(key=lambda r: r['coverage'])
    anchor = None
    if policy_row and policy_row.get('policy_status') == 'target_met':
        n = int(number(policy_row, 'selected_n'))
        require(n in by_n, 'Locked selective-policy point is absent from its saved curve.')
        anchor = by_n[n]
        require(abs(anchor['selective_error']-number(policy_row, 'selective_error')) < 1e-8,
                'Selective-policy error differs from the corresponding saved curve point.')
    candidates = {}
    for percent in range(1, 101):
        target = percent/100
        best = min(possible, key=lambda row: abs(row['coverage']-target))
        if anchor and 0 < abs(best['selected_n']-anchor['selected_n']) < MIN_CURVE_BLOCK:
            continue
        candidates[best['selected_n']] = best
    if anchor:
        candidates[anchor['selected_n']] = anchor
    chosen = []
    for n, row in sorted(candidates.items()):
        if chosen and n-chosen[-1]['selected_n'] < MIN_CURVE_BLOCK:
            continue
        chosen.append(row)
    require(len(chosen) >= 20, 'Too few aggregate selective-curve points remain.')
    require(all(b['selected_n']-a['selected_n'] >= MIN_CURVE_BLOCK for a,b in zip(chosen,chosen[1:])),
            'Selective curve does not meet the retained-count increment rule.')
    if anchor:
        require(anchor in chosen, 'The public curve lost its locked-policy anchor.')
    return chosen

def make_bundle(base_dir=BASE_DIR, notebook_file=NOTEBOOK_FILE, output_parent=OUTPUT_PARENT):
    base = Path(base_dir).expanduser().resolve()
    run = base/'outputs'/'hoi_trust_analysis_v4'
    tables, protocols = run/'tables', run/'protocols'
    manifest = json.loads(source_bytes(protocols/'reproducibility_manifest.json'))
    final = json.loads(source_bytes(protocols/'final_run_record.json'))
    for label, record in [('reproducibility manifest', manifest), ('final-run record', final)]:
        for key, expected in [('analysis_fingerprint', EXPECTED_FINGERPRINT), ('code_version_sha256', EXPECTED_CODE_VERSION),
                              ('direct_package_lock_sha256', EXPECTED_DIRECT_SHA), ('resolved_environment_sha256', EXPECTED_ENV_SHA),
                              ('python_version', EXPECTED_PYTHON)]:
            require(record.get(key) == expected, f'Unexpected {key} in {label}. The saved run is not the reviewed manuscript run.')
    settings = manifest.get('settings', {})
    require(manifest.get('run_profile') == 'definitive' and settings.get('outer_repeats') == 3
            and settings.get('outer_splits') == 5 and settings.get('inner_splits') == 3
            and settings.get('bootstrap') == 2000, 'The completed definitive 3-repeat/2,000-bootstrap run is required.')
    require(final.get('all_cells_reached_final_audit') is True and final.get('submission_artifacts_complete') is True
            and final.get('environment_warnings') == 0, 'The saved final audit did not complete successfully.')
    notebook, direct_lock, notebook_info = snapshot_notebook(base, notebook_file)
    artifacts = {'HOI_TRUST_pipeline.ipynb': notebook}
    provenance = []
    checks = []
    environment_notes = []
    for filename, expected, fallback in [
        ('requirements_direct_lock.txt', EXPECTED_DIRECT_SHA, direct_lock),
        ('requirements_frozen_full.txt', EXPECTED_ENV_SHA,
         ('\n'.join(manifest.get('full_environment', []))+'\n').encode())]:
        path = protocols/filename
        raw = source_bytes(path) if path.is_file() else fallback
        require(sha(raw) == expected, f'{filename} is missing or differs from the recorded environment. '
                'Recover the original file; do not substitute a current pip freeze.')
        require(b'://' not in raw and b' @ ' not in raw, f'Inspect local paths or package URLs before exporting {filename}.')
        artifacts['environment/'+filename] = raw
        environment_notes.append({'file': filename, 'sha256': sha(raw),
                                  'origin': 'saved protocol file' if path.is_file() else 'reconstructed from verified recorded run information'})
    # Only these aggregate fields are exported. Unknown fields are never copied.
    loaded = {}
    def collect(name, required, extra=(), maximum=250, source_dir=tables):
        rows, record = rows_from(source_dir/name, required, maximum, list(required)+list(extra))
        loaded[name] = rows; provenance.append(record)
        artifacts['results/'+name] = csv_bytes(rows)
        return rows

    candidate_extra = ('repeats mean_AUROC sd_AUROC mean_AUPRC sd_AUPRC mean_AUPRC_lift '
        'mean_macro_pathogen_AUROC mean_macro_pathogen_AUPRC_lift mean_precision_at_5pct mean_recall_at_5pct '
        'mean_enrichment_at_5pct scaled_mean_macro_pathogen_AUPRC_lift scaled_mean_recall_at_5pct '
        'scaled_mean_enrichment_at_5pct selection_score').split()
    candidates = collect('candidate_summary.csv', ['pathway', 'candidate'], candidate_extra, 16)
    require(len(candidates) == 16 and {(r['pathway'],r['candidate']) for r in candidates} == {(p,c) for p in PATHWAYS for c in CANDIDATES},
            'Candidate summary must contain eight unique candidates per pathway.')
    require(all(number(r,'repeats') == 3 for r in candidates), 'Candidate results do not all contain three repeats.')
    for p, auroc, auprc in [('patient', .6144, .3887), ('environmental', .6089, .4423)]:
        rounded_check(candidates,'mean_AUROC',auroc,pathway=p,candidate='Strict nested stack')
        rounded_check(candidates,'mean_AUPRC',auprc,pathway=p,candidate='Strict nested stack')
    collect('candidate_metrics_by_repeat.csv', ['pathway','candidate','repeat','n'],
            METRICS+CAPACITY+['macro_pathogen_AUROC','macro_pathogen_AUPRC_lift'], 48)
    repeats = loaded['candidate_metrics_by_repeat.csv']
    require(len(repeats) == 48 and {(r['pathway'],r['candidate'],int(number(r,'repeat'))) for r in repeats}
            == {(p,c,i) for p in PATHWAYS for c in CANDIDATES for i in range(3)}, 'Incomplete candidate/repeat export.')
    absolute = collect('component_bootstrap_95ci.csv',
        ['pathway','model','metric','estimate','ci_95_lower','ci_95_upper','bootstrap_replicates'],
        ['capacity_fraction_for_queue_metrics','resampling_unit'], 20)
    paired = collect('paired_stack_vs_best_single_95ci.csv',
        ['pathway','metric','paired_difference','paired_ci_95_lower','paired_ci_95_upper','bootstrap_replicates'],
        ['comparison','best_single_model','capacity_fraction_for_queue_metrics','stack_estimate','single_estimate','resampling_unit'], 10)
    require(len(absolute) == 20 and len(paired) == 10 and all(number(r,'bootstrap_replicates') == 2000 for r in absolute+paired),
            'Expected complete 2,000-replicate primary and paired summaries.')
    for r in absolute:
        require(number(r,'ci_95_lower') <= number(r,'ci_95_upper'), 'Reversed primary interval.')
    for r in paired:
        require(number(r,'paired_ci_95_lower') <= number(r,'paired_ci_95_upper'), 'Reversed paired interval.')
    for p, auc, ap in [('patient',.6154,.3895),('environmental',.6102,.4436)]:
        rounded_check(absolute,'estimate',auc,pathway=p,model='Strict nested stack',metric='AUROC')
        rounded_check(absolute,'estimate',ap,pathway=p,model='Strict nested stack',metric='AUPRC')
    for p,m,x in [('patient','AUROC',.0023),('patient','AUPRC',.0022),('environmental','AUROC',.0020),('environmental','AUPRC',.0030)]:
        rounded_check(paired,'paired_difference',x,pathway=p,metric=m)
    calibration = collect('crossfitted_internal_calibration_by_repeat.csv', ['pathway','repeat','n','components'], METRICS, 6)
    require({(r['pathway'],int(number(r,'repeat'))) for r in calibration} == {(p,i) for p in PATHWAYS for i in range(3)},
            'Cross-fitted calibration must contain all six pathway/repeat rows.')
    collect('probability_calibrator_apparent_fit.csv', ['pathway'], METRICS+['locked_candidate','calibration_method','partition'], 2)
    policy = collect('v2_locked_operating_policies.csv', ['pathway','policy_type','policy_status'],
        CAPACITY+['selection_rule','confidence_threshold','coverage','selective_error','selective_accuracy','selected_n'], 4)
    for p,precision,recall in [('patient',.4207,.1392),('environmental',.5045,.1489)]:
        rounded_check(policy,'precision_at_capacity',precision,pathway=p,policy_type='capacity')
        rounded_check(policy,'recall_at_capacity',recall,pathway=p,policy_type='capacity')
        require(number(one(policy,pathway=p,policy_type='capacity'),'capacity_fraction') == .10, 'Unexpected locked review capacity.')
    require(one(policy,pathway='environmental',policy_type='selective')['policy_status'] == 'target_not_met',
            'Environmental selective policy must remain infeasible for the manuscript run.')
    sensitivities = collect('strict_sensitivity_results.csv', ['pathway','scenario','model','observations'],
        ['events','dependency_components','capacity_fraction','precision','recall','enrichment',
         'delta_vs_primary_AUROC','delta_vs_primary_AUPRC','delta_vs_primary_precision','delta_vs_primary_recall','delta_vs_primary_enrichment']+METRICS, 18)
    require(sum(r['model']=='Strict nested stack' for r in sensitivities) == 9, 'Expected nine applicable strict-stack sensitivities.')
    collect('pathogen_context_sensitivity.csv', ['pathway','feature_set'],
        ['mean_AUROC','mean_AUPRC','mean_macro_pathogen_AUROC','mean_macro_pathogen_AUPRC_lift','mean_recall_at_5pct'], 4)
    compatibility = collect('v2_legacy_holdout_metrics_NOT_CONFIRMATORY.csv', ['pathway'], METRICS+CAPACITY+['evaluation_status'], 2)
    for p,auc,ap in [('patient',.6113,.3845),('environmental',.6146,.4525)]:
        rounded_check(compatibility,'AUROC',auc,pathway=p);rounded_check(compatibility,'AUPRC',ap,pathway=p)
    collect('v2_legacy_capacity_metrics_NOT_CONFIRMATORY.csv', ['pathway','locked_capacity_fraction','actual_selected_fraction','precision','recall'],
        ['number_needed_to_review'], 2)  # Per-pathogen score-boundary dictionary is deliberately omitted.
    uncertainty = collect('v2_legacy_uncertainty_metrics_NOT_CONFIRMATORY.csv',
        ['pathway','nominal_coverage','observation_coverage','component_complete_coverage','singleton_fraction','ambiguous_fraction','empty_fraction'],
        ['selective_nonabstention_fraction','selective_error_nonabstained'], 2)
    for p,c,a in [('patient',.8962,.8942),('environmental',.9309,.9207)]:
        rounded_check(uncertainty,'component_complete_coverage',c,pathway=p)
        rounded_check(uncertainty,'ambiguous_fraction',a,pathway=p)
    pathogens = collect('legacy_within_pathogen_metrics_NOT_CONFIRMATORY.csv', ['pathway','run','observations','events','nonevents'],
        METRICS+CAPACITY+['evaluation_status','observation_conformal_coverage','component_complete_conformal_coverage'], 22)
    require(len(pathogens) == 22 and all(number(r,'events')>=30 and number(r,'nonevents')>=30 for r in pathogens),
            'Unexpected support or count in pathogen aggregates.')
    collect('v2_locked_candidate_grouped_permutation_importance.csv',
        ['pathway','domain','mean_AUPRC_decrease','lower_2_5','upper_97_5'],
        ['locked_candidate','n_columns','baseline_operating_AUPRC','sd_AUPRC_decrease'], 8)
    collect('v2_stack_meta_coefficients.csv', ['pathway','base_model','meta_logit_coefficient'], maximum=10)
    subgroups, record = rows_from(tables/'v2_legacy_subgroup_robustness_NOT_CONFIRMATORY.csv',
        ['pathway','subgroup','value','n','cases','controls'], 100,
        ['pathway','subgroup','value','n','cases','controls','adequate_support']+METRICS+CAPACITY+['conformal_coverage','abstention_fraction'])
    subgroups = [r for r in subgroups if number(r,'cases')>=30 and number(r,'controls')>=30]
    require(bool(subgroups), 'No eligible aggregate subgroup rows remain.')
    record['exported_rows'] = len(subgroups); provenance.append(record)
    artifacts['results/v2_legacy_subgroup_robustness_NOT_CONFIRMATORY.csv'] = csv_bytes(subgroups)
    collect('v2_partition_governance.csv', ['pathway','partition','observations','dependency_components'], maximum=10)
    partitions = collect('five_way_partition_summary_blinded.csv', ['pathway','partition','observations','cases','dependency_components'],
        ['observations','cases','controls','case_proportion','patients','matched_sets','dependency_components'], 10, base/'outputs'/'tables')
    counts = []
    for row in partitions:
        row = dict(row)
        if row['partition'] == 'final_test':
            row['cases'] = sum(int(number(r,'events')) for r in pathogens if r['pathway'] == row['pathway'])
        row['cases'] = int(number(row,'cases'))
        row['controls'] = int(number(row,'observations'))-row['cases']
        row['case_proportion'] = row['cases']/number(row,'observations')
        row['partition'] = {'final_test':'compatibility','operating_point':'operating_policy'}.get(row['partition'],row['partition'])
        counts.append(row)
    for pathway, expected_n, expected_cases in [('patient',56305,16933),('environmental',43403,14923)]:
        require(sum(int(number(r,'observations')) for r in counts if r['pathway']==pathway)==expected_n
                and sum(r['cases'] for r in counts if r['pathway']==pathway)==expected_cases, 'Analytical-role totals differ from the frozen cohort.')
    artifacts['results/analytical_partition_counts.csv'] = csv_bytes(counts)
    collect('pathway_identity_summary.csv', ['analysis','observations','cases','pathway_specific_patients','matched_sets','dependency_components'],
        ['irregular_matched_sets'], 2, base/'outputs'/'tables')
    collect('pathogen_pathway_summary.csv', ['analysis','run','observations','cases','controls','case_proportion'],
        ['pathway_specific_patients'], 22, base/'outputs'/'tables')
    for p in sorted(PATHWAYS):
        collect(f'capacity_curve_operating_{p}.csv', ['capacity_fraction','recall_at_capacity'], CAPACITY, 5)
        private, record = rows_from(tables/f'selective_risk_coverage_operating_{p}.csv',
            ['coverage','selective_error','selected_n'], 100_000,
            ['coverage','selective_error','selected_n'])
        public = coarsen_curve(private, one(policy,pathway=p,policy_type='selective'))
        record.update(exported_rows=len(public), transformation='Coarse coverage grid, no confidence scores; minimum retained-count increment 30')
        provenance.append(record)
        artifacts[f'results/selective_curve_public_{p}.csv'] = csv_bytes(public)
    dashboard, record = rows_from(tables/'dashboard_aggregate_summary.csv',
        ['pathway','pathogen_code','uncertainty','observations','selected_for_review'], 100,
        ['pathway','pathogen_code','uncertainty','observations','selected_for_review'])
    queue, comp = {}, {}
    for row in dashboard:
        p, code, category = row['pathway'],row['pathogen_code'],row['uncertainty']
        require(code in PATHOGENS and category in {'Singleton','Ambiguous','Empty'}, 'Unexpected dashboard category.')
        n, selected = int(number(row,'observations')),int(number(row,'selected_for_review'))
        require(0 <= selected <= n, 'Invalid aggregate dashboard review count.')
        queue.setdefault((p,code), [0,0]);queue[p,code][0]+=n;queue[p,code][1]+=selected
        comp[p,category] = comp.get((p,category),0)+n
    require(sum(v[0] for v in queue.values()) == 19850 and sum(v[1] for v in queue.values()) == 1995,
            'Dashboard totals differ from the manuscript (19,850 scored; 1,995 selected).')
    artifacts['results/dashboard_review_queue.csv'] = csv_bytes([
        {'pathway':p,'pathogen_code':c,'observations':v[0],'selected_for_review':v[1]} for (p,c),v in sorted(queue.items())])
    artifacts['results/dashboard_uncertainty_composition.csv'] = csv_bytes([
        {'pathway':p,'uncertainty':c,'observations':n} for (p,c),n in sorted(comp.items())])
    record['transformation'] = 'Exported two aggregate margins; original score means and three-way cells omitted'
    provenance.append(record)
    # The following embedded script reads only the public aggregate CSVs.
    artifacts['scripts/regenerate_aggregate_figures.py'] = FIGURE_SCRIPT.encode()
    artifacts['environment/python_version.txt'] = (EXPECTED_PYTHON+'\n').encode()
    self_path = globals().get('__file__')
    if self_path and Path(self_path).is_file():
        self_raw = Path(self_path).read_bytes()
        if EXPORTER_MARKER.encode() in self_raw:
            artifacts['scripts/export_hoi_release.py'] = self_raw
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    public_manifest = {
        'project':'HOI-TRUST','release_version':RELEASE_VERSION,'exported_utc':stamp,
        'analysis_fingerprint':EXPECTED_FINGERPRINT,'code_version_sha256':EXPECTED_CODE_VERSION,
        'run_profile':manifest['run_profile'],'settings':settings,'python_version':EXPECTED_PYTHON,
        'direct_package_lock_sha256':EXPECTED_DIRECT_SHA,'resolved_environment_sha256':EXPECTED_ENV_SHA,
        'original_execution_completed_utc':final.get('execution_completed_utc'),
        'notebook':notebook_info,'environment':environment_notes,'source_aggregates':provenance,
        'model_fitting_performed_by_exporter':False,
        'checks':['Reviewed run identity and completed audit','Exact 52-code-cell identity',
                  'Recorded environment hashes','Manuscript primary numerical anchors',
                  'Explicit aggregate column selection','Subgroup minimum support',
                  'Selective-curve count increments','Dashboard aggregate totals'],
        'limitations':['Source CSVs are checked against the saved run and numerical anchors; this is not a new statistical validation.',
                       'Fine-grained selective confidence thresholds are not public; Figure 3 uses the documented coarse aggregate grid.',
                       'Source data, prediction rows, fitted models and caches are excluded.'],
    }
    artifacts['release_manifest.json'] = (json.dumps(public_manifest,indent=2,ensure_ascii=False)+'\n').encode()
    artifacts['README.md'] = README_TEXT.encode()
    artifacts['LICENSE'] = LICENSE_TEXT.encode()
    artifacts['SHA256SUMS.txt'] = ''.join(f'{sha(raw)}  {name}\n' for name,raw in sorted(artifacts.items())).encode()
    parent = Path(output_parent).expanduser().resolve() if output_parent else base.parent/'HOI_TRUST_public_release'
    require(parent != base and run not in parent.parents and parent != run, 'Release output must be outside the analysis-run folder.')
    parent.mkdir(parents=True,exist_ok=True)
    name = f'HOI_TRUST_{RELEASE_VERSION}_release_candidate_{stamp}'
    destination, archive = parent/name, parent/(name+'.zip')
    require(not destination.exists() and not archive.exists(), 'Release output already exists.')
    staging = Path(tempfile.mkdtemp(prefix='hoi_export_',dir=parent))
    try:
        for relative, raw in artifacts.items():
            target = staging/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        for relative, raw in artifacts.items():
            require(sha((staging/relative).read_bytes()) == sha(raw), f'Export verification failed: {relative}')
        staging.rename(destination)
        with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for relative in sorted(artifacts):
                z.write(destination/relative,arcname=relative)
        with zipfile.ZipFile(archive) as z:
            require(z.testzip() is None and set(z.namelist()) == set(artifacts), 'ZIP validation failed.')
            for relative, raw in artifacts.items():
                require(sha(z.read(relative)) == sha(raw), f'ZIP content mismatch: {relative}')
    except Exception:
        if staging.exists():shutil.rmtree(staging)
        if archive.exists():archive.unlink()
        if destination.exists():shutil.rmtree(destination)
        raise
    print('\nEXPORT CHECKS PASSED')
    print(f'Files: {len(artifacts)}; aggregate CSVs: {sum(n.endswith(".csv") for n in artifacts)}')
    print(f'Release folder: {destination}')
    print(f'Release ZIP: {archive}')
    print('Unzip and inspect before uploading. Upload the contents to GitHub, then create a new release.')
    print('Regenerate figures with: python scripts/regenerate_aggregate_figures.py')
    return {'directory':str(destination),'zip':str(archive),'files':len(artifacts)}

def main():
    if str(BASE_DIR).startswith('/content/drive/') and not Path('/content/drive/MyDrive').is_dir():
        try:
            from google.colab import drive
            drive.mount('/content/drive')
        except ImportError:
            raise ExportError('Mount Google Drive in Colab, or set BASE_DIR to your local saved project folder.')
    result = make_bundle()
    if DOWNLOAD_ZIP:
        try:
            from google.colab import files
            files.download(result['zip'])
        except ImportError:
            pass
    return result

# Embedded aggregate plotting script and release documentation are inserted below.

FIGURE_SCRIPT = '"""Rebuild aggregate HOI-TRUST figures without patient-level data or model fitting.\n\nFigure 3 uses a documented coarse grid with a minimum retained-count increment.\nThese are aggregate-only reproductions; fine-grained private curve points are omitted.\n"""\nfrom pathlib import Path\nimport argparse\nimport numpy as np\nimport pandas as pd\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib.ticker import PercentFormatter, MaxNLocator\n\nCOLORS = {\'patient\': \'#df7900\', \'environmental\': \'#2462a5\'}\nORDER = [\'patient\', \'environmental\']\nPATHOGENS = {\n    \'DS_E_coli\': \'Drug-susceptible E. coli\', \'ESBL_E_coli\': \'ESBL E. coli\',\n    \'DS_K_pneumoniae\': \'Drug-susceptible K. pneumoniae\',\n    \'ESBL_K_pneumoniae\': \'ESBL K. pneumoniae\', \'C_diff\': \'C. difficile\',\n    \'VSE_faecalis\': \'Vancomycin-susceptible E. faecalis\',\n    \'VRE_faecium\': \'Vancomycin-resistant E. faecium\', \'MSSA\': \'MSSA\', \'MRSA\': \'MRSA\',\n    \'DS_P_aeruginosa\': \'Drug-susceptible P. aeruginosa\',\n    \'DR_P_aeruginosa\': \'Drug-resistant P. aeruginosa\',\n}\n\ndef render_all(root):\n    root = Path(root)\n    tables, figures = root / \'results\', root / \'figures\'\n    figures.mkdir(exist_ok=True)\n    plt.rcParams.update({\'font.size\': 11, \'axes.spines.top\': False,\n                         \'axes.spines.right\': False, \'savefig.dpi\': 300})\n\n    def read(name):\n        return pd.read_csv(tables / name)\n\n    def save(fig, name):\n        fig.savefig(figures / (name + \'.png\'), bbox_inches=\'tight\', facecolor=\'white\')\n        fig.savefig(figures / (name + \'.pdf\'), bbox_inches=\'tight\', facecolor=\'white\')\n        plt.close(fig)\n\n    # Primary component-bootstrap estimates and paired differences.\n    absolute = read(\'component_bootstrap_95ci.csv\')\n    paired = read(\'paired_stack_vs_best_single_95ci.csv\')\n    fig, ax = plt.subplots(1, 3, figsize=(16, 6), layout=\'constrained\',\n                           gridspec_kw={\'width_ratios\': [1, 1, 1.3]})\n    combinations = [(\'patient\', \'Strict nested stack\'), (\'patient\', \'Best single model: CatBoost\'),\n                    (\'environmental\', \'Strict nested stack\'), (\'environmental\', \'Best single model: Extra Trees\')]\n    labels = [\'Patient — stack\', \'Patient — CatBoost\', \'Environmental — stack\', \'Environmental — Extra Trees\']\n    for j, metric in enumerate([\'AUROC\', \'AUPRC\']):\n        for y, (pathway, model) in enumerate(combinations):\n            row = absolute.loc[(absolute.pathway == pathway) & (absolute.model == model)\n                               & (absolute.metric == metric)].iloc[0]\n            x, lo, hi = row[\'estimate\'], row.ci_95_lower, row.ci_95_upper\n            ax[j].errorbar(x, y, xerr=[[x-lo], [hi-x]], fmt=\'o\' if \'stack\' in model else \'s\',\n                           color=COLORS[pathway], capsize=3)\n            ax[j].annotate(f\'{x:.4f}\', (x, y), xytext=(6, -12), textcoords=\'offset points\', fontsize=9)\n        ax[j].set(title=metric, xlim=(0, 1), yticks=range(4))\n        ax[j].set_yticklabels(labels if j == 0 else [])\n        ax[j].invert_yaxis()\n        ax[j].set_xlabel(\'Estimate (95% component-bootstrap interval)\')\n        ax[j].grid(axis=\'x\', alpha=.15)\n    metrics = [\'AUROC\', \'AUPRC\', \'precision\', \'recall\', \'enrichment\']\n    for pathway, shift in [(\'patient\', -.11), (\'environmental\', .11)]:\n        for y, metric in enumerate(metrics):\n            row = paired.loc[(paired.pathway == pathway) & (paired.metric == metric)].iloc[0]\n            x, lo, hi = row.paired_difference, row.paired_ci_95_lower, row.paired_ci_95_upper\n            ax[2].errorbar(x, y+shift, xerr=[[x-lo], [hi-x]], fmt=\'o\', color=COLORS[pathway],\n                           capsize=3, label=pathway.title() if y == 0 else None)\n    ax[2].axvline(0, color=\'.4\', ls=\'--\')\n    ax[2].set(title=\'Paired differences\', yticks=range(5),\n              yticklabels=[\'AUROC\', \'AUPRC\', \'Precision at 10%\', \'Recall at 10%\', \'Enrichment at 10%\'],\n              xlabel=\'Stack minus pathway-specific comparator\')\n    ax[2].invert_yaxis(); ax[2].legend(frameon=False)\n    fig.suptitle(\'Primary internal validation: estimates and paired differences\')\n    save(fig, \'results_primary_component_bootstrap\')\n\n    # Single-repeat sensitivities; uncertainty intervals were not estimated.\n    data = read(\'strict_sensitivity_results.csv\')\n    data = data.loc[data.model == \'Strict nested stack\']\n    scenarios = [\'one_exact_copy\', \'adults_only\', \'exclude_irregular_matched_sets\',\n                 \'exclude_resolved_any_abx_summary\', \'exclude_largest_dependency_component\']\n    names = [\'One exact copy\', \'Adults only\', \'Exclude irregular matched sets\',\n             \'Remove antibiotic summary (patient only)\', \'Exclude largest component\']\n    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5), sharey=True, layout=\'constrained\')\n    for axis, key, title in zip(axes, [\'delta_vs_primary_AUROC\', \'delta_vs_primary_AUPRC\',\n                                     \'delta_vs_primary_precision\'], [\'Change in AUROC\', \'Change in AUPRC\',\n                                                                     \'Change in precision at 10%\']):\n        for pathway, shift in [(\'patient\', -.10), (\'environmental\', .10)]:\n            sub = data.loc[data.pathway == pathway].set_index(\'scenario\')\n            for y, scenario in enumerate(scenarios):\n                if scenario in sub.index:\n                    axis.scatter(sub.loc[scenario, key], y+shift, color=COLORS[pathway],\n                                 label=pathway.title() if y == 0 else None)\n        axis.axvline(0, color=\'.4\', ls=\'--\'); axis.set_title(title)\n        axis.set_xlabel(\'Difference from corresponding primary repeat\')\n        axis.xaxis.set_major_locator(MaxNLocator(4))\n        axis.grid(axis=\'x\', alpha=.15)\n    axes[0].set(yticks=range(5), yticklabels=names); axes[0].invert_yaxis()\n    axes[-1].legend(frameon=False)\n    fig.suptitle(\'Development sensitivities: single-repeat diagnostics without confidence intervals\')\n    save(fig, \'results_sensitivity_analyses\')\n\n    # Capacity points and coarsened selective error curve, with no per-record thresholds.\n    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout=\'constrained\')\n    for pathway in ORDER:\n        cap = read(f\'capacity_curve_operating_{pathway}.csv\').sort_values(\'capacity_fraction\')\n        curve = read(f\'selective_curve_public_{pathway}.csv\').sort_values(\'coverage\')\n        axes[0].plot(cap.capacity_fraction, cap.recall_at_capacity, \'-o\', color=COLORS[pathway], label=pathway.title())\n        axes[1].plot(curve.coverage, curve.selective_error, color=COLORS[pathway], label=pathway.title())\n    axes[0].axhline(.10, color=\'.4\', ls=\'--\'); axes[0].axvline(.10, color=\'.4\', ls=\':\')\n    axes[1].axhline(.20, color=\'.4\', ls=\'--\')\n    for axis in axes:\n        axis.xaxis.set_major_formatter(PercentFormatter(1))\n        axis.yaxis.set_major_formatter(PercentFormatter(1)); axis.legend(frameon=False)\n        axis.grid(alpha=.15)\n    axes[0].set(title=\'Case capture by review capacity\', xlabel=\'Nominal within-pathogen review fraction\', ylabel=\'Case recall\')\n    axes[1].set(title=\'Selective error on the public aggregate grid\', xlabel=\'Retained fraction\', ylabel=\'Classification error\')\n    fig.suptitle(\'Operating-policy development; the 20% capacity point is descriptive\')\n    save(fig, \'operating_policy_curves\')\n\n    # Component and observation coverage and set composition.\n    uncertainty = read(\'v2_legacy_uncertainty_metrics_NOT_CONFIRMATORY.csv\').set_index(\'pathway\').loc[ORDER]\n    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout=\'constrained\')\n    x = np.arange(2)\n    for shift, key, label, color in [(-.16, \'observation_coverage\', \'Observation\', \'#4c78a8\'),\n                                     (.16, \'component_complete_coverage\', \'Complete component\', \'#df7900\')]:\n        bars = axes[0].bar(x+shift, uncertainty[key], .30, label=label, color=color)\n        axes[0].bar_label(bars, labels=[f\'{v:.4f}\' for v in uncertainty[key]], padding=3)\n    axes[0].axhline(.90, ls=\'--\', color=\'.4\', label=\'Nominal 90%\')\n    bottom = np.zeros(2)\n    for key, label, color in [(\'singleton_fraction\', \'Singleton\', \'#4c78a8\'),\n                               (\'ambiguous_fraction\', \'Ambiguous\', \'#f2a43c\'), (\'empty_fraction\', \'Empty\', \'#e45756\')]:\n        v = uncertainty[key].to_numpy()\n        axes[1].bar(x, v, bottom=bottom, label=label, color=color)\n        for j in range(2):\n            if v[j] > .001: axes[1].text(j, bottom[j]+v[j]/2, f\'{v[j]:.2%}\', ha=\'center\', va=\'center\')\n        bottom += v\n    for axis, title in zip(axes, [\'Empirical coverage\', \'Prediction-set composition\']):\n        axis.set(xticks=x, xticklabels=[p.title() for p in ORDER], title=title)\n        axis.yaxis.set_major_formatter(PercentFormatter(1))\n        axis.legend(frameon=False, loc=\'upper center\', bbox_to_anchor=(.5, -.10), ncol=3, fontsize=9)\n    axes[0].set_ylim(0, 1.09); axes[1].set_ylim(0, 1)\n    fig.suptitle(\'Conformal uncertainty in the non-confirmatory compatibility partition\')\n    save(fig, \'results_conformal_compatibility\')\n\n    # Pathogen-specific lift.\n    pathogens = read(\'legacy_within_pathogen_metrics_NOT_CONFIRMATORY.csv\')\n    fig, axis = plt.subplots(figsize=(12, 8), layout=\'constrained\')\n    for pathway, shift in [(\'patient\', -.10), (\'environmental\', .10)]:\n        sub = pathogens.loc[pathogens.pathway == pathway].set_index(\'run\')\n        for y, code in enumerate(PATHOGENS):\n            row = sub.loc[code]\n            axis.scatter(row.AUPRC_lift, y+shift, s=max(20, float(row.observations)/12),\n                         color=COLORS[pathway], alpha=.85, label=pathway.title() if y == 0 else None)\n    axis.axvline(1, color=\'.4\', ls=\'--\')\n    axis.set(yticks=range(11), yticklabels=list(PATHOGENS.values()),\n             xlabel=\'AUPRC / pathogen-specific sampled case proportion\',\n             title=\'Within-pathogen AUPRC lift: non-confirmatory compatibility partition\')\n    axis.invert_yaxis(); axis.legend(frameon=False)\n    save(fig, \'results_within_pathogen_heterogeneity\')\n\n    # Grouped permutation percentiles are not confidence intervals.\n    reliance = read(\'v2_locked_candidate_grouped_permutation_importance.csv\')\n    domains = [\'Colonisation pressure\', \'Pathogen context and interactions\', \'Comorbidity\',\n               \'Demographic and surgery\', \'Prior antibiotic exposure\']\n    fig, axis = plt.subplots(figsize=(12, 6), layout=\'constrained\')\n    for pathway, shift in [(\'patient\', -.12), (\'environmental\', .12)]:\n        sub = reliance.loc[reliance.pathway == pathway].set_index(\'domain\')\n        labelled = False\n        for y, domain in enumerate(domains):\n            if domain not in sub.index: continue\n            row = sub.loc[domain]; x = row.mean_AUPRC_decrease\n            axis.errorbar(x, y+shift, xerr=[[x-row.lower_2_5], [row.upper_97_5-x]], fmt=\'o\',\n                          color=COLORS[pathway], capsize=4, label=pathway.title() if not labelled else None)\n            labelled = True\n    axis.axvline(0, color=\'.4\', ls=\'--\')\n    axis.set(yticks=range(5), yticklabels=domains, xlabel=\'AUPRC decrease after joint domain permutation\',\n             title=\'Grouped permutation reliance on the operating partition\')\n    axis.invert_yaxis(); axis.legend(frameon=False)\n    fig.supxlabel(\'Bars are empirical permutation percentiles, not confidence intervals. Associations are not causal.\', fontsize=9)\n    save(fig, \'results_grouped_permutation_reliance\')\n\n    # Dashboard figure uses two aggregate margins, without individual scores.\n    queue = read(\'dashboard_review_queue.csv\')\n    composition = read(\'dashboard_uncertainty_composition.csv\')\n    fig, axes = plt.subplots(1, 2, figsize=(14, 7), layout=\'constrained\')\n    for pathway, shift in [(\'patient\', -.16), (\'environmental\', .16)]:\n        sub = queue.loc[queue.pathway == pathway].set_index(\'pathogen_code\')\n        vals = [sub.loc[code, \'selected_for_review\'] for code in PATHOGENS]\n        axes[0].barh(np.arange(11)+shift, vals, .30, color=COLORS[pathway], label=pathway.title())\n    axes[0].set(yticks=range(11), yticklabels=list(PATHOGENS.values()), xlabel=\'Selected observations\', title=\'Review queue by pathogen target\')\n    axes[0].invert_yaxis(); axes[0].legend(frameon=False)\n    bottom = np.zeros(2)\n    totals = composition.groupby(\'pathway\').observations.sum()\n    for label, color in [(\'Singleton\', \'#4c78a8\'), (\'Ambiguous\', \'#f2a43c\'), (\'Empty\', \'#e45756\')]:\n        sub = composition.loc[composition.uncertainty == label].set_index(\'pathway\')\n        vals = np.array([sub.loc[p, \'observations\']/totals[p] if p in sub.index else 0 for p in ORDER])\n        axes[1].bar(np.arange(2), vals, bottom=bottom, label=label, color=color); bottom += vals\n    axes[1].set(xticks=range(2), xticklabels=[p.title() for p in ORDER], title=\'Uncertainty composition\', ylim=(0, 1))\n    axes[1].yaxis.set_major_formatter(PercentFormatter(1)); axes[1].legend(frameon=False)\n    selected = int(queue.selected_for_review.sum()); n = int(queue.observations.sum())\n    ambiguous = int(composition.loc[composition.uncertainty == \'Ambiguous\', \'observations\'].sum())\n    fig.suptitle(f\'Retrospective surveillance review — {n:,} observations; {selected:,} selected ({selected/n:.1%}); {ambiguous/n:.1%} ambiguous\')\n    save(fig, \'results_dashboard_demonstrator\')\n    print(f\'Saved seven aggregate figure reproductions to {figures}\')\n\nif __name__ == \'__main__\':\n    parser = argparse.ArgumentParser(description=__doc__)\n    parser.add_argument(\'--root\', type=Path, default=Path(__file__).resolve().parents[1])\n    args = parser.parse_args()\n    render_all(args.root)\n'

README_TEXT = "# HOI-TRUST\n\nExplainable and uncertainty-aware machine learning for hospital-onset infection surveillance.\n\nThis release contains analysis source, recorded environment specifications and aggregate results\nfor the completed retrospective internal-validation run. It provides a research prototype;\nexternal, genuinely temporal and prospective clinical validation have not been performed.\nThe compatibility partition is non-confirmatory because it informed earlier development.\n\n## Contents\n\n- `HOI_TRUST_pipeline.ipynb`: the original 52 code-cell sources, with saved outputs, cell metadata,\n  execution counts, attachments and widget state removed.\n- `results/`: aggregate numerical tables. Counts and metrics are conditional on the matched sample.\n- `environment/`: the recorded direct dependency lock, full environment freeze and Python version.\n- `scripts/regenerate_aggregate_figures.py`: plotting from the public aggregate CSVs only.\n- `scripts/export_hoi_release.py`: the export script, when executed from its standalone `.py` file.\n- `release_manifest.json`: run identity, source-export hashes, transformations and checks.\n- `SHA256SUMS.txt`: checksums for the other release files.\n\nThe manifest's code-cell-source hash uses newline-joined cell sources. It is distinct from the\noriginal notebook's embedded code-version hash, which uses the pipeline's own hashing convention.\nThe stripped notebook also has a different whole-file hash from the retained executed notebook.\n\n## Source data and analysis execution\n\nThe source dataset is PhysioNet's *Predictors of Hospital Onset Infection*, version 1.0.0:\nhttps://doi.org/10.13026/k70x-0m81. Credentialing, required training and the applicable data-use\nagreement are required. Source records, individual predictions, models and private caches are\nnot included in this release.\n\nTo run the original analysis, use Python **3.13.16** and the recorded dependencies. Obtain authorised\nsource files and configure `HOI_TRUST_DATA_DIR` for a private project folder. Save a copy of the\nnotebook there under its expected name, `HOI_TRUST_submission_ready.ipynb`, and set\n`HOI_TRUST_NOTEBOOK_PATH` to that copy before executing the configuration cell. Follow the run\ninstructions in the notebook, including the package-installation restart. Use\n`FAST_MODE=False` and `HOI_TRUST_PROFILE=definitive` for the reported analysis. The notebook\nchecks source-release counts and data, code, environment and cache identities. Source-specific\nfingerprints can differ on a different authorised release. Do not bypass those checks to mix runs.\n\nThe exporter does not require the current interpreter to match the modelling environment: it\ncopies or reconstructs the **recorded** environment and verifies its original hashes. It never\nuses a new `pip freeze` as a substitute for the historical environment.\n\n## Regenerate the public aggregate figures\n\nFrom the extracted release root, with NumPy, pandas and Matplotlib available:\n\n```bash\npython scripts/regenerate_aggregate_figures.py\n```\n\nSeven PNG/PDF aggregate reproductions will be written to `figures/`. This command reads no\nindividual records, loads no fitted models and performs no model fitting or performance resampling.\nThese plotting scripts were added for the release; the original analysis notebook remains available\nfor the original plotting procedure. Layouts may differ from the submitted figure files.\n\n| Manuscript item | Aggregate numerical sources in `results/` |\n|---|---|\n| Table 1 | `candidate_summary.csv`, `candidate_metrics_by_repeat.csv` |\n| Table 2 and Figure 1 | `component_bootstrap_95ci.csv`, `paired_stack_vs_best_single_95ci.csv` |\n| Table 3 | `probability_calibrator_apparent_fit.csv`, `v2_locked_operating_policies.csv` |\n| Figure 2 | `strict_sensitivity_results.csv` (strict nested stack rows only) |\n| Figure 3 | `capacity_curve_operating_patient.csv`, `capacity_curve_operating_environmental.csv`, `selective_curve_public_patient.csv`, `selective_curve_public_environmental.csv` |\n| Figure 4 | `v2_legacy_uncertainty_metrics_NOT_CONFIRMATORY.csv` |\n| Figure 5 | `legacy_within_pathogen_metrics_NOT_CONFIRMATORY.csv` |\n| Figure 6 | `v2_locked_candidate_grouped_permutation_importance.csv` |\n| Figure 7 | `dashboard_review_queue.csv`, `dashboard_uncertainty_composition.csv` |\n| Main cohort/partition tables | `pathway_identity_summary.csv`, `analytical_partition_counts.csv` |\n| Appendix calibration/subgroups/context sensitivity | `crossfitted_internal_calibration_by_repeat.csv`, `v2_legacy_subgroup_robustness_NOT_CONFIRMATORY.csv`, `pathogen_context_sensitivity.csv` |\n\nThe saved compatibility metrics CSV contains fresh 5%-capacity diagnostics; transferred locked\n10%-capacity results are in `v2_legacy_capacity_metrics_NOT_CONFIRMATORY.csv`. Do not substitute\none for the other. Repeat-specific calibration estimates and candidate repeat means can differ\nfrom the primary estimate formed by averaging repeated predictions by observation before resampling.\n\n## Public aggregation and interpretation\n\nThe exporter selects named aggregate fields only. Small subgroup rows with fewer than 30 cases\nor 30 controls are omitted. Dashboard data are exported as two margins: pathway/pathogen review\ncounts and pathway/uncertainty counts. Individual scores, score means and the original three-way\ndashboard cells are omitted.\n\nThe selective-error curve is shared on an approximately 1%-spaced coverage grid with at least\n30 retained observations and a minimum difference of 30 retained observations between adjacent\npublic points. Feasible locked-policy points are preserved. Per-point confidence thresholds are\nnot exported. Values are selected from the saved aggregate curve without interpolation; the\nvery small retained-fraction tail and fine-grained geometry of the submitted curve are omitted.\nExact fine-grained reconstruction requires authorised source access and the analysis notebook.\nThe 30-observation rule is an export convention, not a clinical performance guarantee or formal\nproof of anonymisation.\n\nSensitivity estimates are single-repeat diagnostics without dedicated confidence intervals.\nPermutation bars describe perturbation percentiles, not confidence intervals. Subgroups and\ncompatibility results are descriptive. Matched-sample probabilities are not population absolute risks.\n\nThe exporter verifies the saved run identifiers, the source notebook code, recorded environment\nhashes, selected manuscript numerical anchors, CSV structure and public output checksums. This\ndoes not independently establish the provenance of every CSV or reproduce the statistical analysis.\nInspect the release candidate before publication. Preserve the original restricted run separately.\n\n## Version and citation\n\nThis package is prepared for release **v1.0.1**; the notebook's internal analysis version remains\n3.1.0. After publishing and archiving the release, add its version-specific Zenodo DOI and correct\nsoftware-author metadata here and use that DOI in the manuscript. A release-candidate ZIP does\nnot establish that the files are already public.\n\n## Licence\n\nThe author-generated code is MIT licensed, as in the existing repository. This does not replace\nthe source dataset's access conditions or grant permission to redistribute restricted records.\n"

LICENSE_TEXT = 'MIT License\n\nCopyright (c) 2026 Joel-AWC\n\nPermission is hereby granted, free of charge, to any person obtaining a copy\nof this software and associated documentation files (the "Software"), to deal\nin the Software without restriction, including without limitation the rights\nto use, copy, modify, merge, publish, distribute, sublicense, and/or sell\ncopies of the Software, and to permit persons to whom the Software is\nfurnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\nOUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\nSOFTWARE.\n'

if __name__ == "__main__":
    try:
        main()
    except ExportError as exc:
        print("EXPORT STOPPED:", exc)
        raise

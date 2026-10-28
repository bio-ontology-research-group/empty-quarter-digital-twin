"""Independent pathway audit, matched ecological contrasts and numerical sensitivity."""
from dataclasses import replace
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

from eq_biodt.v03.evidence import source_directory
from eq_biodt.v03.process import Parameters, State, simulate
from eq_biodt.v03.planning import parameter_ensemble, compare
from eq_biodt.v03.interventions import intervene

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from independent_reference import module_definition, satisfied


def cluster_ci(frame, columns, statistic, rng, draws=2000):
    """Resample sites with replacement, preserving all rows within each site."""
    sums = frame.groupby('site')[columns].sum().to_numpy()
    counts = frame.groupby('site').size().to_numpy()
    idx = rng.integers(0, len(sums), (draws, len(sums)))
    means = sums[idx].sum(axis=1) / counts[idx].sum(axis=1)[:,None]
    estimates = statistic(means)
    return [float(x) for x in np.quantile(estimates,[.025,.975])]


def analyse(root, output, backtest, scenarios):
    rng = np.random.default_rng(42)
    mag = pd.read_csv(root/'data/v03/mag_evidence.tsv', sep='\t')
    panel = pd.read_csv(root/'data/v03/sample_evidence.tsv', sep='\t')
    expr = module_definition((source_directory()/'M00165.txt').read_text())
    raw = pd.read_csv(source_directory()/'rep_ko_v2.tsv.gz',sep='\t',names=['rep','kos'],keep_default_na=False)
    calls = raw.set_index('rep').kos.map(lambda s: satisfied(expr,set(s.split(','))))
    actual = mag.set_index('rep').cbb_all_steps_detected
    disagreement = int((calls.reindex(actual.index)!=actual).sum())
    if disagreement:
        raise ValueError('Independent module evaluation disagrees with production evidence')
    ontology = json.loads((root/'data/v03/entailed_memberships.json').read_text())
    owl_calls = pd.Series({k: 'CBBCompleteMAG' in v for k,v in ontology['classes_by_mag'].items()})
    owl_disagreement = int((calls.reindex(owl_calls.index)!=owl_calls).sum())
    if owl_disagreement:
        raise ValueError('Independent module evaluation disagrees with OWL entailments')
    columns = [c for c in panel if 'cbb' in c]
    contrasts=[]
    for endpoint in ['cbb_all_steps_detected_share_mapped','h2_cbb_all_steps_share_mapped']:
        if endpoint not in panel:
            raise ValueError(f'Unknown endpoint {endpoint}; available: {columns}')
        grouped=panel.groupby(['site','trip_num','compartment'])[endpoint].mean().unstack()
        for a,b in [('Deep','PlantRoot'),('Deep','Surface'),('PlantRoot','Surface')]:
            paired=grouped[[a,b]].dropna().reset_index()
            paired['difference']=100*(paired[a]-paired[b])
            contrasts.append({'endpoint':endpoint,'contrast':f'{a} minus {b}',
                'matched_site_campaigns':len(paired),'sites':int(paired.site.nunique()),
                'mean_difference_percentage_points':float(paired.difference.mean()),
                'site_bootstrap_95_percentile':cluster_ci(paired,['difference'],lambda x:x[:,0],rng)})
    pred=pd.DataFrame(backtest['predictions'])
    pred['squared_error']=(pred.predicted_log1p_ge-pred.observed_log1p_ge)**2
    wide=pred.pivot(index=['site','compartment','campaign'],columns='model',values='squared_error').reset_index()
    comparisons=[]
    for baseline in ['persistence','compartment_mean','global_mean']:
        cols=['geography_season_ridge',baseline]
        comparisons.append({'baseline':baseline,'visits':len(wide),'sites':int(wide.site.nunique()),
          'rmse_difference_ridge_minus_baseline':float(np.sqrt(wide[cols[0]].mean())-np.sqrt(wide[cols[1]].mean())),
          'site_bootstrap_95_percentile':cluster_ci(wide,cols,lambda x:np.sqrt(x[:,0])-np.sqrt(x[:,1]),rng)})
    coverage=[]
    for (model,compartment),g in pred.groupby(['model','compartment']):
        coverage.append({'model':model,'compartment':compartment,'visits':len(g),'coverage':float(g.covered_90.mean())})
    p=Parameters(k_h2_s=.001,k_co_s=.001)
    convergence=[]
    for scenario in scenarios:
        values={str(n):simulate(State(1.8),scenario['forcing'],replace(p,integration_steps_per_day=n))['summary']['h2_uptake_mol_m2'] for n in (24,96,384)}
        convergence.append({'scenario':scenario['name'],'h2_uptake_by_steps':values,
          'relative_error_96_vs_384':abs(values['96']-values['384'])/values['384'] if values['384'] else 0.})
    fine = compare(State(1.8), scenarios, parameter_ensemble(replace(p, integration_steps_per_day=384)))
    coarse = compare(State(1.8), scenarios, parameter_ensemble(p))
    rank_stable = [r['rank_reversal'] for r in fine['pairwise_uptake_comparisons']] == [r['rank_reversal'] for r in coarse['pairwise_uptake_comparisons']]
    if not rank_stable:
        raise ValueError('Scenario rank reversal is unresolved at default integration resolution')
    grid_relative_errors = [abs(a['h2_uptake_mol_m2']-b['h2_uptake_mol_m2'])/b['h2_uptake_mol_m2'] for a,b in zip(coarse['outcomes'],fine['outcomes']) if b['h2_uptake_mol_m2']]
    # Same-input flat model check: wrapper must change only the specified leaf.
    forcing=scenarios[0]['forcing']
    query={'forcing':{'0':{'water_mm':12}}}
    paired=intervene(State(1.8),forcing,p,query,root/'data/v03')
    flat=[dict(f) for f in forcing];flat[0]['water_mm']=12
    assert paired['intervened']==simulate(State(1.8),flat,p)
    result={'seed':42,'bootstrap_draws':2000,'pathway_audit':{'mags':len(mag),'full_module_candidates':int(calls.sum()),
         'flat_vs_production_disagreements':disagreement,'flat_vs_owl_disagreements':owl_disagreement,
         'transaldolase_route_without_sbpase':int((actual & ~mag.set_index('rep').sbpase_detected).sum()),
         'archaeal_rubisco_with_prk':int((mag.domain.eq('Archaea') & mag.rbcl_prk_detected).sum()),
         'scope':'Independent parser of pinned KEGG expression, same KO information; no experimental truth labels.'},
       'matched_compartment_contrasts':contrasts,'forecast_paired_comparisons':comparisons,
       'coverage_by_compartment':coverage,'numerical_convergence':convergence,
       'intervention_flat_equivalence':True,
       'grid_rank_reversals_stable_96_vs_384':rank_stable,
       'grid_max_relative_integration_error':max(grid_relative_errors),
       'uncertainty_scope':'Site resampling retains all within-site visits; intervals describe this observational cohort. No treatment randomization or adjustment for site selection.'}
    (output/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result

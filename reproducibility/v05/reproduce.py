"""Reproduce the v0.5 manuscript from bundled sources (offline after dependency and reasoner installation)."""
from pathlib import Path
import subprocess,sys,json,hashlib
ROOT=Path(__file__).resolve().parents[2]
def execute(relative,*args):subprocess.run([sys.executable,str(ROOT/relative),*args],cwd=ROOT,check=True)
def main():
 for record in json.loads((ROOT/'data/external/pmc-supplement-manifest.json').read_text()):
  path=ROOT/'data/external'/record['study']/record['file']
  body=path.read_bytes()
  if len(body)!=record['bytes'] or hashlib.sha256(body).hexdigest()!=record['sha256'] or hashlib.md5(body).hexdigest()!=record['md5']:
   raise ValueError(f"External source digest mismatch: {path}")
 for key,info in json.loads((ROOT/'data/sources/manifest.json').read_text())['sources'].items():
  path=ROOT/'data/sources'/info['file']
  if hashlib.sha256(path.read_bytes()).hexdigest()!=info['sha256']:
   raise ValueError(f"Source digest mismatch: {key}")
 execute('reproducibility/v03/reproduce.py','--paper','paper')
 execute('reproducibility/v05/external_transfer.py')
 execute('reproducibility/v05/primary_connection.py')
 execute('reproducibility/v05/audit.py')
 execute('reproducibility/v05/figures.py')
 out=ROOT/'reviews/v05';tables=ROOT/'paper/tables'
 out.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
 r=json.loads((out/'external_transfer.json').read_text())
 labels={'frozen_reference':'Frozen reference','calibrated_kinetics':'Kinetic coefficient calibrated','calibrated_diffusion_only':'Diffusion-only coefficient calibrated','training_mean':'Training mean','linear_moisture_temperature':'Linear moisture and temperature','calibrated_retention':'Soil-response parameters calibrated'}
 metrics={row['model']:row for row in r['sonoran']['field_metrics']}
 (tables/'v05_external.tex').write_text('\n'.join(f"{label} & {metrics[key]['rmse']:.3f} & {metrics[key]['mae']:.3f} \\\\" for key,label in labels.items())+'\n')
 import pandas as pd
 rows=pd.read_csv(out/'primary_campaigns.tsv',sep='\t')
 (tables/'v05_campaigns.tex').write_text('\n'.join(f"{row.trip_num} & {row.samples} & {row.sites} & {row.first_date} & {row.last_date} \\\\" for row in rows.itertuples())+'\n')
 import numpy as np
 raw=pd.read_csv(out/'cross_desert_marker_comparison.tsv',sep='\t')
 # Report each fraction's median, never sum component medians into a composition.
 med=raw.groupby('environment')[['h2_1h_fraction','h2_1l_fraction','h2_2a_fraction']].median()*100
 (tables/'v05_markers.tex').write_text('\n'.join(f"{name} & {len(raw[raw.environment.eq(name)])} & {row.iloc[0]:.1f} & {row.iloc[1]:.1f} & {row.iloc[2]:.1f} \\\\" for name,row in med.iterrows())+'\n')
 inputs=[ROOT/'data/external'/x['study']/x['file'] for x in json.loads((ROOT/'data/external/pmc-supplement-manifest.json').read_text())]
 inputs += [ROOT/'reviews/v05/external-evaluation-protocol.md',ROOT/'data/sources/manifest.json']
 inputs += sorted((ROOT/'reproducibility/v05').glob('*.py'))
 digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 outputs=sorted(out.glob('*.tsv'))+[out/'deduction_causal_trace.json',out/'external_transfer.json',out/'independent_audit.json',out/'primary_connection.json']+sorted(tables.glob('v05_*.tex'))
 manifest={'release':'0.5.0','inputs':{str(p.relative_to(ROOT)):digest(p) for p in inputs},'analysis_outputs':{str(p.relative_to(ROOT)):digest(p) for p in outputs},'scope':'Offline analysis from released primary annotations and bundled external supplements after dependency and reasoner installation; raw-read processing is upstream.'}
 (out/'reproduction_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('v0.5 reproduction complete; independent numerical audit passed')
if __name__=='__main__':main()

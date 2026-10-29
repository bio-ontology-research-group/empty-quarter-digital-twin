"""Independent numerical and denominator checks for the external analyses.

Does not call the production gas-exchange function. Recomputes held-out fluxes
from the published equation and verifies native source counts and aggregation.
"""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.optimize import minimize_scalar
from itertools import product
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'reviews/v05'

def check(cond,msg="Audit check failed"):
 if not cond:
  raise ValueError(msg)

def reference(theta,temp,k,phi=.4,b=4.,entry=-.001,q10=2.,diffusion_only=False):
 theta=np.asarray(theta,float);temp=np.asarray(temp,float)
 psi=entry*(theta/phi)**(-b)
 a=np.clip(np.log(85/abs(psi))/np.log(850),0,1)
 if diffusion_only:a=np.ones_like(a)
 d=7.7e-5*(phi-theta)**(10/3)/phi**2
 reaction=k*a*q10**((temp-25)/10)
 return 2.13e-5*np.sqrt(d*reaction)*np.tanh(.1*np.sqrt(reaction/d))*1e9

def run():
 for record in json.loads((ROOT/'data/external/pmc-supplement-manifest.json').read_text()):
  path=ROOT/'data/external'/record['study']/record['file']
  body=path.read_bytes()
  check(len(body)==record['bytes'] and hashlib.sha256(body).hexdigest()==record['sha256'] and hashlib.md5(body).hexdigest()==record['md5'],f"External source digest mismatch: {path}")
 for key,info in json.loads((ROOT/'data/sources/manifest.json').read_text())['sources'].items():
  path=ROOT/'data/sources'/info['file']
  check(hashlib.sha256(path.read_bytes()).hexdigest()==info['sha256'],f"Source digest mismatch: {key}")
 report=json.loads((OUT/'external_transfer.json').read_text())
 obs=pd.read_csv(OUT/'sonoran_field_observations.tsv',sep='\t')
 pred=pd.read_csv(OUT/'sonoran_field_predictions.tsv',sep='\t')
 raw=pd.read_csv(ROOT/'data/external/buzzard2022/peerj-10-12966-s008.csv')
 raw['campaign']=pd.to_datetime(raw.Date).dt.strftime('%Y-%m-%d')
 group=-raw.groupby(['campaign','Treatment'])['Hydrogen.Uptake'].mean()
 for row in obs.itertuples():check(np.isclose(row.y,group.loc[(row.campaign,row.Treatment)]),f"Sonoran mean mismatch: {row.campaign} {row.Treatment}")
 check(len(pred)==16*6 and not pred.duplicated(['campaign','plot','model']).any(),"Invalid Sonoran prediction rows")
 discrepancies=[]
 all_campaigns=set(obs.campaign.unique())
 for fit in report['sonoran']['fits']:
  target=fit['heldout_campaign']
  check(target not in fit['training_campaigns'] and set(fit['training_campaigns'])==all_campaigns-{target},"Sonoran campaign holdout boundary violated")
  test=obs[obs.campaign.eq(target)];train=obs[obs.campaign.ne(target)]
  opt_kin=minimize_scalar(lambda lk:np.mean((reference(train.theta,train.temp,10**lk)-train.y)**2),bounds=(-9,2),method='bounded',options={'xatol':1e-9})
  opt_dif=minimize_scalar(lambda lk:np.mean((reference(train.theta,train.temp,10**lk,diffusion_only=True)-train.y)**2),bounds=(-9,2),method='bounded',options={'xatol':1e-9})
  check(np.isclose(10**opt_kin.x,fit['kinetics_k'],rtol=1e-6,atol=1e-9) and np.isclose(10**opt_dif.x,fit['diffusion_k'],rtol=1e-6,atol=1e-9),"Independent single-parameter fit mismatch")
  expected={
   'frozen_reference':reference(test.theta,test.temp,.001),
   'calibrated_kinetics':reference(test.theta,test.temp,fit['kinetics_k']),
   'calibrated_diffusion_only':reference(test.theta,test.temp,fit['diffusion_k'],diffusion_only=True),
   'training_mean':np.repeat(train.y.mean(),len(test)),
   'linear_moisture_temperature':np.column_stack([np.ones(len(test)),test.theta,test.temp])@np.linalg.lstsq(np.column_stack([np.ones(len(train)),train.theta,train.temp]),train.y,rcond=None)[0]}
  p=fit['tuned'];expected['calibrated_retention']=reference(test.theta,test.temp,p['k'],p['phi'],p['b'],p['entry'],p['q10'])
  # Refit the grid with an independent vectorized equation using training rows only.
  best=float('inf')
  for phi,b,entry,q10 in product([.35,.4,.5],[2.,4.,6.],[-.0001,-.001,-.01],[1.,2.]):
   loss=lambda lk:np.mean((reference(train.theta,train.temp,10**lk,phi,b,entry,q10)-train.y)**2)
   optimum=minimize_scalar(loss,bounds=(-9,2),method='bounded')
   best=min(best,float(optimum.fun))
  check(np.isclose(best,fit['tuned_training_mse'],rtol=1e-6,atol=1e-9),"Independent grid training MSE mismatch")
  scored=np.mean((reference(train.theta,train.temp,p['k'],p['phi'],p['b'],p['entry'],p['q10'])-train.y)**2)
  check(np.isclose(scored,best,rtol=1e-6,atol=1e-9),"Selected grid member training MSE mismatch")
  for model,values in expected.items():
   saved=pred[pred.campaign.eq(target)&pred.model.eq(model)].set_index('plot').loc[test.Treatment,'predicted'].to_numpy()
   discrepancy=float(np.max(abs(values-saved)));discrepancies.append(discrepancy)
   check(np.allclose(values,saved,rtol=1e-10,atol=1e-10),f"Prediction mismatch: {target} {model} {discrepancy}")
 for metric in report['sonoran']['field_metrics']:
  g=pred[pred.model.eq(metric['model'])]
  check(np.isclose(metric['rmse'],np.sqrt(np.mean((g.observed-g.predicted)**2))) and np.isclose(metric['mae'],np.mean(np.abs(g.observed-g.predicted))),f"Field metric mismatch: {metric['model']}")
 lab=pd.read_csv(OUT/'sonoran_lab_observations.tsv',sep='\t')
 check(len(lab)==48 and lab.groupby('source_batch').size().eq(3).all(),"Sonoran lab source-batch count mismatch")
 namib=pd.read_csv(OUT/'namib_observations.tsv',sep='\t')
 check(len(namib)==60 and namib.groupby(['site','state']).size().eq(5).all(),"Namib observation count mismatch")
 check(namib[namib.state.eq('dry')].site.nunique()==2,"Namib quantified dry site count mismatch")
 check(np.allclose(namib.k_per_hour*.0047436*530e-9*1e9/2,namib.nmol_h_g,rtol=5e-6,atol=0),"Namib unit conversion mismatch")
 check(all(c['n_dry_rates']==0 and 'dry_mean' not in c and 'fold' not in c for c in report['namib']['contrasts'][2:]),"Unquantified dry Namib sites must not have imputed dry_mean or fold")
 namib_contrasts=pd.read_csv(OUT/'namib_hydration_contrasts.tsv',sep='\t')
 check(namib_contrasts.loc[namib_contrasts.n_dry_rates.eq(0),['dry_mean','fold','fold_lcl','fold_ucl']].isna().all().all(),"Unquantified dry Namib TSV rows must remain NA")
 temp=pd.read_csv(OUT/'temperate_leave_soil_out.tsv',sep='\t')
 check(len(temp)==77*3 and not temp.duplicated(['heldout_soil','sample_id','replicate','model']).any(),"Temperate prediction row count mismatch")
 # Intentionally outside field fitting: a distinct endpoint in native bottle units.
 rates=pd.read_csv(ROOT/'reviews/v03/external_bottle_rates.tsv',sep='\t')
 for soil in rates.soil.unique():
  train=rates[rates.soil.ne(soil)];test=rates[rates.soil.eq(soil)]
  for model in ['diffusion_only','activation_diffusion','training_mean']:
   if model=='training_mean':values=np.repeat(train.bottle_k_min.mean(),len(test))
   else:
    x=train[model].to_numpy();coef=max(0,float(np.dot(x,train.bottle_k_min)/np.dot(x,x)));values=test[model].to_numpy()*coef
   saved=temp[temp.heldout_soil.eq(soil)&temp.model.eq(model)]
   check(list(zip(test.sample_id,test.replicate))==list(zip(saved.sample_id,saved.replicate)),"Temperate replicate ordering mismatch")
   check(np.allclose(values,saved.predicted_k_min,rtol=1e-10,atol=1e-10),"Temperate leave-soil-out prediction mismatch")
 trace=json.loads((OUT/'deduction_causal_trace.json').read_text())
 check(trace['ontology_supported']['fluxes']==trace['same_information_flat']['fluxes'],"Same-information flat resource trace mismatch")
 overlap=pd.read_csv(OUT/'capacity_overlap_bounds.tsv',sep='\t').dropna()
 check(len(overlap)==550 and ((overlap.genome_resolved>=overlap.lower-1e-12)&(overlap.genome_resolved<=overlap.upper+1e-12)).all(),"Frechet overlap bounds violated")
 result={'status':'passed' ,'independently_recomputed_field_predictions':len(pred),'max_absolute_field_prediction_discrepancy':max(discrepancies),
         'sonoran_lab_bottles':len(lab),'namib_quantified_rates':len(namib),'temperate_prediction_rows':len(temp),
         'checks':['Corrected native source endpoint and sign','All source-batch conditions retained','Whole-campaign holdout boundaries and independent training-only grid refit','Independent diffusion-reaction equation','All field metrics','Unquantified dry Namib endpoints retained','Complete temperate denominator']}
 (OUT/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':run()

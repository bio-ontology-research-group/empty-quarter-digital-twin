"""Trace primary desert samples through deductive capacity and causal admission."""
from pathlib import Path
from contextlib import closing
import json,hashlib,warnings
import numpy as np,pandas as pd,openpyxl
from eq_biodt.v03.api import resource
from eq_biodt.v03.resources import resource_ceiling
from eq_biodt.v03.ontology import load_classification
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'reviews/v05'

def verify_digest(path,expected_sha256,expected_md5=None,expected_bytes=None):
 body=path.read_bytes()
 if hashlib.sha256(body).hexdigest()!=expected_sha256:
  raise ValueError(f"SHA-256 mismatch: {path}")
 if expected_md5 is not None and hashlib.md5(body).hexdigest()!=expected_md5:
  raise ValueError(f"MD5 mismatch: {path}")
 if expected_bytes is not None and len(body)!=expected_bytes:
  raise ValueError(f"Byte-length mismatch: {path}")
 return path

def run():
 OUT.mkdir(parents=True,exist_ok=True)
 v03_manifest=json.loads((ROOT/'data/v03/manifest.json').read_text())['artifacts']
 for name,digest in v03_manifest.items():
  verify_digest(ROOT/'data/v03'/name,digest)
 src_manifest=json.loads((ROOT/'data/sources/manifest.json').read_text())['sources']
 verify_digest(ROOT/'data/sources'/src_manifest['primary_field_ledger']['file'],src_manifest['primary_field_ledger']['sha256'])
 ext_records=json.loads((ROOT/'data/external/pmc-supplement-manifest.json').read_text())
 jordaan_rec=next(r for r in ext_records if r['study']=='jordaan2020' and r['file']=='mSystems.01131-20-st004.xlsx')
 verify_digest(ROOT/'data/external/jordaan2020/mSystems.01131-20-st004.xlsx',jordaan_rec['sha256'],jordaan_rec['md5'],jordaan_rec['bytes'])
 panel=pd.read_csv(ROOT/'data/v03/sample_evidence.tsv',sep='\t')
 mags=pd.read_csv(ROOT/'data/v03/mag_evidence.tsv',sep='\t')
 visits=panel.groupby(['trip_num','site'],as_index=False).size()
 sites=panel.groupby('site',as_index=False).agg(latitude=('lat','first'),longitude=('lon','first'),samples=('sample','size'),campaigns=('trip_num','nunique'))
 sites.to_csv(OUT/'primary_sites.tsv',sep='\t',index=False)
 field=pd.read_csv(ROOT/'data/sources/trip1-field-ledger.tsv',sep='\t')
 field['site']=pd.to_numeric(field.site,errors='coerce');field=field.dropna(subset=['site'])
 xy=field.coordinates.str.extract(r'([0-9.]+) N,\s*([0-9.]+) E').astype(float)
 field['field_latitude']=xy[0];field['field_longitude']=xy[1]
 joined=sites[sites.site.le(60)].merge(field[['site','field_latitude','field_longitude']],on='site',validate='one_to_one')
 if len(joined)!=59 or not (np.allclose(joined.latitude,joined.field_latitude,rtol=0,atol=1e-8) and np.allclose(joined.longitude,joined.field_longitude,rtol=0,atol=1e-8)):
  raise ValueError("Core site coordinates do not match field ledger")
 joined.to_csv(OUT/'primary_field_coordinate_audit.tsv',sep='\t',index=False)
 campaigns=panel.groupby('trip_num',as_index=False).agg(samples=('sample','size'),sites=('site','nunique'),first_date=('collection_date','min'),last_date=('collection_date','max'))
 campaigns.to_csv(OUT/'primary_campaigns.tsv',sep='\t',index=False)
 # Sharp Frechet overlap bounds from community marginals, same mapped universe.
 h=panel.h2_carrier_union_share_mapped;c=panel.cbb_all_steps_detected_share_mapped
 overlap=pd.DataFrame({'sample':panel['sample'],'site':panel.site,'h2_marginal':h,'cbb_marginal':c,
                       'lower':np.maximum(0,h+c-1),'upper':np.minimum(h,c),'genome_resolved':panel.h2_cbb_all_steps_share_mapped})
 valid=overlap.dropna()
 if not ((valid.genome_resolved>=valid.lower-1e-12)&(valid.genome_resolved<=valid.upper+1e-12)).all():
  raise ValueError("Genome-resolved intersection violates Frechet bounds")
 overlap.to_csv(OUT/'capacity_overlap_bounds.tsv',sep='\t',index=False)
 # Fixed illustrative supplies: a worked conditional query, not field observations.
 values=json.loads((ROOT/'data/v03/example_resource.json').read_text())
 chosen=mags.loc[mags.h2_cbb_all_steps & ~mags.co_form1].sort_values('rep').iloc[0]
 values['mag_id']=chosen.rep
 exact=resource(values)
 supplies={k:v for k,v in values.items() if k!='mag_id'}
 flat=resource_ceiling(**supplies,pathway_supported=bool(chosen.h2_cbb_all_steps),co_oxidation_supported=bool(chosen.co_form1))
 conflated=resource_ceiling(**supplies,pathway_supported=True,co_oxidation_supported=True)
 if exact['fluxes']!=flat['fluxes'] or exact['net_species_exchange_mol_m2_s']!=flat['net_species_exchange_mol_m2_s']:
  raise ValueError("Ontology-admitted resource ceiling disagrees with same-information flat calculation")
 used=load_classification(ROOT/'data/v03')
 trace={'inputs':values,'declared_scenario':'Illustrative supplied substrate/energy limits; no field-flux calibration',
        'mag_annotations':{k:chosen[k].item() if hasattr(chosen[k],'item') else chosen[k] for k in ['rep','cbb_kos','h2_1l','h2_1h','h2_2a','co_form1','completeness','contamination']},
        'entailed_classes':used['classes_by_mag'][chosen.rep],'ontology_supported':exact,'same_information_flat':flat,'unsupported_CO_admission':conflated,
        'scope':'The contrast tests information admission. Missing CO entailment excludes that reaction from the declared model, without establishing biological absence.'}
 (OUT/'deduction_causal_trace.json').write_text(json.dumps(trace,indent=2,allow_nan=False)+'\n')
 report={'samples':len(panel),'usable_mapping_samples':int(panel.mapped_read_pct.notna().sum()),'sites':len(sites),
         'core_sites':int(sites.site.le(60).sum()),'auxiliary_sites':sites.loc[sites.site.gt(60),'site'].tolist(),
         'site_campaigns':len(visits),'compartments':panel.compartment.value_counts().to_dict(),
         'latitude_range':[float(panel.lat.min()),float(panel.lat.max())],'longitude_range':[float(panel.lon.min()),float(panel.lon.max())],
         'mapping_overlap_samples':len(valid),'ambiguous_marginals':int((valid.upper-valid.lower>1e-12).sum()),
         'median_overlap_interval_percentage_points':float(np.median(valid.upper-valid.lower)*100),
         'co_admission_ceiling_ratio':conflated['fluxes']['carbohydrate_c']/exact['fluxes']['carbohydrate_c'],
         'data_boundary':'Primary samples and derived annotations/read summaries. No whole-desert spatial interpolation or raw-read reassembly is claimed.',
         'campaign_dates':'Curated metadata snapshot; fourth-campaign dates are 2024-08-29 to 2024-08-31, differing from a companion draft narrative.'}
 (OUT/'primary_connection.json').write_text(json.dumps(report,indent=2)+'\n')
 # Descriptive external genomic comparison: same three marker families, no union claim.
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',category=UserWarning,module=r'openpyxl\.worksheet\._reader')
  with closing(openpyxl.load_workbook(ROOT/'data/external/jordaan2020/mSystems.01131-20-st004.xlsx',read_only=True,data_only=True)) as w:
   vals=list(w['Annotation summary'].values)
 header=list(vals[1]);rows=[]
 for row in vals[2:14]:
  if row[0] not in ['C1','Namib','Gobi','Mojave']:continue
  marker=np.array([row[header.index(k)] for k in ['Hyd-1h','Hyd-1l','Hyd-2a']],float)/100
  rows.append({'environment':{'C1':'South Australian desert'}.get(row[0],row[0]),'sample':row[0],**dict(zip(['h2_1h_ge','h2_1l_ge','h2_2a_ge'],marker)),
               'normalizer':'Jordaan single-copy ribosomal marker normalization','source':'10.1128/mSystems.01131-20 Table S4'})
 for r in panel.itertuples():
  marker=np.array([r.hydB_1h,r.hydB_1l,r.hydB_2a],float)
  if not np.isfinite(marker).all():continue
  rows.append({'environment':"Rub al-Khali",'sample':r.sample,**dict(zip(['h2_1h_ge','h2_1l_ge','h2_2a_ge'],marker)),
               'normalizer':'Rub al-Khali screened single-copy marker normalization','source':'Bundled primary read-screen table'})
 comparison=pd.DataFrame(rows)
 cols=['h2_1h_ge','h2_1l_ge','h2_2a_ge']
 for col in cols:comparison[col.replace('_ge','_fraction')]=comparison[col]/comparison[cols].sum(axis=1).replace(0,np.nan)
 comparison.to_csv(OUT/'cross_desert_marker_comparison.tsv',sep='\t',index=False)
 print(json.dumps(report,indent=2))
if __name__=='__main__':run()

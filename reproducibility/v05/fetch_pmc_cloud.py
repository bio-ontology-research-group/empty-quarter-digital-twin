"""Retrieve supplementary data through PMC's 2026 public cloud interface."""
from pathlib import Path
import concurrent.futures as cf,httpx,hashlib,json
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'data/external'
STUDIES={'tribbia2026':'PMC13488325','jordaan2020':'PMC7677003','buzzard2022':'PMC8934528'}
def fetch(item):
 study,pmc=item;base='https://pmc-oa-opendata.s3.amazonaws.com/'
 url=f'{base}{pmc}.1/{pmc}.1.json'
 r=httpx.get(url,timeout=30);r.raise_for_status();meta=r.json()
 (OUT/study/'pmc-cloud.json').write_bytes(r.content)
 print(study,'metadata',str(meta)[:1000],flush=True)
 return study,meta
with cf.ThreadPoolExecutor(max_workers=3) as ex:
 results=list(ex.map(fetch,STUDIES.items()))
from urllib.parse import urlparse,parse_qs
jobs=[]
for study,meta in results:
 for uri in meta['media_urls']:
  parsed=urlparse(uri);name=Path(parsed.path).name
  if not any(name.lower().endswith(ext) for ext in ('.xlsx','.csv','.rmd','.r','.docx','.pdf')):continue
  jobs.append((study,name,uri,meta['license_code']))
def download(job):
 study,name,uri,license=job;parsed=urlparse(uri)
 url='https://pmc-oa-opendata.s3.amazonaws.com'+parsed.path
 r=httpx.get(url,timeout=30);r.raise_for_status();body=r.content
 expected=parse_qs(parsed.query)['md5'][0]
 assert hashlib.md5(body).hexdigest()==expected,(study,name,'MD5 mismatch')
 (OUT/study/name).write_bytes(body)
 return {'study':study,'file':name,'url':url,'license':license,'sha256':hashlib.sha256(body).hexdigest(),'md5':expected,'bytes':len(body)}
with cf.ThreadPoolExecutor(max_workers=5) as ex:
 records=list(ex.map(download,jobs))
(OUT/'pmc-supplement-manifest.json').write_text(json.dumps(records,indent=2)+'\n')
for r in records:print(r['study'],r['file'],r['bytes'])

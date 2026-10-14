"""Explicitly acquire only confirmed source files from a supplied scope manifest.

No invocation occurred during finalization. List-only is the default. Calling
--download is an explicit user action; it never grants project derivative rights.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,csv,json,urllib.request,urllib.parse,hashlib
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--config',type=Path,required=True)
p.add_argument('--output-root',type=Path,required=True)
p.add_argument('--study',choices=['Berg_ALI','Qin_cotton_rat','HLCA'])
p.add_argument('--download',action='store_true')
a=p.parse_args();rows=[]
with a.config.open() as f:records=list(csv.DictReader(f,delimiter='\t'))
for r in records:
    if a.study and r['study_id']!=a.study:continue
    if r['license']!='CC BY 4.0':raise ValueError('Source-specific acquisition scope not confirmed.')
    u=urllib.parse.urlparse(r['url'])
    if u.scheme!='https' or u.hostname not in ['zenodo.org','datasets.cellxgene.cziscience.com']:raise ValueError('Unexpected source URL.')
    rel=Path(r['source_relpath'])
    if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe relative source path.')
    dest=a.output_root.resolve()/rel
    row=dict(study_id=r['study_id'],url=r['url'],file=str(dest),status='LIST_ONLY_NO_DOWNLOAD',derivative_permission='NOT_CONFERRED')
    if a.download:
        if dest.exists():raise FileExistsError('Refusing to overwrite: '+str(dest))
        dest.parent.mkdir(parents=True,exist_ok=True);tmp=dest.with_suffix(dest.suffix+'.partial')
        if tmp.exists():raise FileExistsError('Prior partial download exists: '+str(tmp))
        md5=hashlib.md5();sha=hashlib.sha256();n=0
        with urllib.request.urlopen(r['url']) as f,tmp.open('xb') as out:
            for b in iter(lambda:f.read(8*1024**2),b''):out.write(b);md5.update(b);sha.update(b);n+=len(b)
        expected=r['official_checksum']
        if expected not in ['','NULL'] and expected!='md5:'+md5.hexdigest():raise ValueError('Official source checksum mismatch; partial file retained.')
        if r.get('bytes') not in [None,'','NULL'] and n!=int(float(r['bytes'])):raise ValueError('Official source byte count mismatch; partial file retained.')
        tmp.rename(dest);row.update(status='DOWNLOADED_VERIFIED_SOURCE_FILE',bytes=n,sha256=sha.hexdigest())
    rows.append(row)
print(json.dumps(rows,ensure_ascii=False,indent=2))

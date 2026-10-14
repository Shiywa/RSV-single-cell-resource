"""Portable source H5AD and unified metadata access with explicit release root."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import json,argparse,sqlite3
import numpy as np,pandas as pd,pyarrow.parquet as pq,anndata as ad

def metadata_root(root):
    root=Path(root)
    return root/'curated_metadata' if (root/'curated_metadata/metadata/cell_index.parquet').is_file() else root

def cells(root,source=None,columns=None):
    p=metadata_root(root)/'metadata/cell_index.parquet';names=pq.read_schema(p).names
    sourcefield='source_label' if 'source_label' in names else 'study_id'
    return pq.read_table(p,columns=columns,filters=[(sourcefield,'=',source)] if source else None).to_pandas()

def marker_counts(root,source,cell_uids,gene_ids):
    root=Path(root);m=cells(root,source)
    selected=m[m.cell_uid.isin(cell_uids)].sort_values('expression_row_index')
    if set(selected.cell_uid)!=set(cell_uids):raise ValueError('Missing requested Cell UID')
    files=selected.source_h5ad_file.unique();assert len(files)==1
    a=ad.read_h5ad(root/files[0],backed='r')
    rows=selected.expression_row_index.to_numpy(np.int64);x=a.X[rows,:].tocsr()
    assert np.array_equal(a.obs_names[rows].to_numpy(),selected.cell_uid)
    assert np.array_equal(a.obs.original_cell_id.iloc[rows].astype(str).to_numpy(),selected.original_cell_id)
    genes=a.var.native_gene_id.astype(str).to_numpy();gi=pd.Index(genes).get_indexer(gene_ids)
    with np.load(root/a.uns['measurement_mask_reference']) as mask:
        assert np.array_equal(mask['genes'],genes)
        byqc=dict(zip(mask['libraries'],mask['availability']));out=[]
        for pos,r in enumerate(selected.itertuples()):
            for gene,j in zip(gene_ids,gi):
                measured=bool(j>=0 and byqc[r.qc_group_id][j])
                out.append(dict(cell_uid=r.cell_uid,original_cell_id=r.original_cell_id,source_id=r.source_id,source_label=source,native_gene_id=gene,measured=measured,count=int(x[pos,j]) if measured else None,source_h5ad_file=r.source_h5ad_file,expression_row_index=int(r.expression_row_index),qc_group_id=r.qc_group_id))
    a.file.close();return pd.DataFrame(out)

def comparisons(root,focused=False):
    p=metadata_root(root)/'comparison/comparison_eligibility.tsv'
    if not p.is_file():return dict(status='BLOCKED_BY_RELEASE_SCOPE',reason='Full project comparison/eligibility assets not included in source-limited candidate.')
    f=pd.read_csv(p,sep='\t',na_values=['NULL'],keep_default_na=False)
    if focused:f=f[f.gate1_5_audited.astype(str).str.lower().eq('true')]
    return f

def identity(root,source):
    p=metadata_root(root)/'metadata/resource.sqlite'
    if not p.is_file():return dict(status='BLOCKED_BY_RELEASE_SCOPE',reason='Complete project identity/evidence database is not included; original source-author labels remain separate.')
    c=sqlite3.connect('file:'+str(p.resolve())+'?mode=ro',uri=True)
    f=pd.read_sql_query('select study_id,source_id,gate2b_canonical_identity,gate2b_fine_subtype,gate2b_identity_status,gate2b_review_status,gate2b_version,count(*) cells_n from annotations_gate2b where study_id=? group by 1,2,3,4,5,6,7',c,params=[source]);c.close();return f

def trace(root,cell_uid):
    p=metadata_root(root)/'metadata/resource.sqlite'
    if p.is_file():
        c=sqlite3.connect('file:'+str(p.resolve())+'?mode=ro',uri=True);c.row_factory=sqlite3.Row
        out=[dict(r) for r in c.execute('select * from publication_expression_mapping where cell_uid=?',[cell_uid])];c.close();return out
    m=cells(root);return m[m.cell_uid.eq(cell_uid)].to_dict('records')

def main():
    p=argparse.ArgumentParser();p.add_argument('--release-root',required=True,type=Path);p.add_argument('--query',choices=['cells','identity','comparisons','trace'],required=True);p.add_argument('--source');p.add_argument('--cell');p.add_argument('--focused',action='store_true');p.add_argument('--limit',type=int,default=10);a=p.parse_args()
    if a.query=='cells':out=cells(a.release_root,a.source).head(a.limit)
    elif a.query=='identity':out=identity(a.release_root,a.source)
    elif a.query=='comparisons':out=comparisons(a.release_root,a.focused)
    else:out=trace(a.release_root,a.cell)
    if isinstance(out,pd.DataFrame):out=out.head(a.limit).where(pd.notna(out),None).to_dict('records')
    print(json.dumps(out,default=str,ensure_ascii=False,indent=2))
if __name__=='__main__':main()

"""Export an authorised certified native count snapshot to a source-local H5AD.

For original-repository reconstruction first use rebuild_original_sources.py
with legitimately available frozen selection/maps and source files. This wrapper
is the same tested source-local projection/streaming writer used in packaging.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd,pyarrow.parquet as pq
from h5ad_io import write_stream

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['native-h5ad','cell-metadata','source-var','measurement-mask','output']:p.add_argument('--'+arg,required=True,type=Path)
    p.add_argument('--source-label',required=True);p.add_argument('--count-technology',required=True);p.add_argument('--mask-reference',required=True)
    p.add_argument('--compression',choices=['gzip','lzf'],default='gzip');p.add_argument('--gzip-level',type=int,default=4);a=p.parse_args()
    m=pq.read_table(a.cell_metadata).to_pandas()
    if 'study_id' in m:m=m[m.study_id.eq(a.source_label)]
    fields=['cell_uid','original_cell_id','source_id','study_id','species_namespace','qc_group_id','matrix_id','native_row_index','sample_id','library_id','verified_unit_id','unit_kind']
    obs=m[[f for f in fields if f in m]].sort_values('native_row_index').copy()
    obs['count_technology']=a.count_technology;obs.index=pd.Index(obs.cell_uid,name='cell_key')
    var=pq.read_table(a.source_var).to_pandas().set_index('feature_key')
    if 'supplementary_homology_symbol' in var:var=var.drop(columns=['supplementary_homology_symbol'])
    columns=var.original_feature_index.to_numpy(np.int64);sp=str(obs.species_namespace.iloc[0])
    uns=dict(resource_version='1.0.0-publication-packaging-rc1',source_label=a.source_label,source_id=sorted(obs.source_id.unique()),expression_type='authorised certified native counts; '+a.count_technology,feature_namespace=sp,annotation_version_if_included='NOT_INCLUDED',provenance_reference='supplied authorised source inventory',measurement_mask_reference=a.mask_reference)
    result=write_stream(a.native_h5ad,a.output,obs,var,columns,a.measurement_mask,uns,a.compression,a.gzip_level)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()

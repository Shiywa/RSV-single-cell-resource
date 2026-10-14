"""Independent source-H5AD/metadata tests. Audit hook precedes third-party imports."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,json,hashlib,os,time,platform,csv
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--release-root',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--arrow-path',type=Path)
p.add_argument('--historical-resource-root',type=Path);p.add_argument('--allowed-work-root',type=Path)
p.add_argument('--forbidden-root',action='append',type=Path,default=[])
p.add_argument('--relocated',action='store_true');a=p.parse_args()
root=a.release_root.resolve();output=a.output.resolve();output.mkdir(parents=True,exist_ok=True)
denied=[];reads=set()
def audit(event,values):
    if event not in ['open','sqlite3.connect','os.listdir','os.scandir'] or not values:return
    x=values[0]
    if not isinstance(x,(str,bytes,os.PathLike)):return
    s=os.fsdecode(x)
    if s.startswith('file:'):s=s[5:].split('?',1)[0]
    q=Path(s).absolute()
    restricted=any(q.is_relative_to(z.resolve()) for z in a.forbidden_root)
    if a.historical_resource_root and q.is_relative_to(a.historical_resource_root.resolve()):
        restricted |= not (a.allowed_work_root and q.is_relative_to(a.allowed_work_root.resolve()))
    if restricted:denied.append(dict(event=event,path=str(q)));raise PermissionError('Forbidden data-path access: '+str(q))
    if event in ['open','sqlite3.connect']:reads.add(str(q))
sys.addaudithook(audit)
if a.arrow_path:sys.path.insert(0,str(a.arrow_path.resolve()))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np,pandas as pd,h5py,anndata as ad,pyarrow as pa,pyarrow.parquet as pq
from load_resource import cells,marker_counts,comparisons,identity,trace,metadata_root
begin=time.perf_counter();checks=[];hashes={}
public=not (root/'curated_metadata').is_dir();scope='PUBLIC_SOURCE_LIMITED' if public else 'INTERNAL_COMPLETE'
env=json.dumps(dict(python=platform.python_version(),anndata=ad.__version__,h5py=h5py.__version__,numpy=np.__version__,pandas=pd.__version__,pyarrow=pa.__version__,host=platform.node(),physical_independent_lab=False))
def digest(file):
    file=Path(file)
    if str(file) not in hashes:
        h=hashlib.sha256()
        with file.open('rb') as f:
            for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
        hashes[str(file)]=h.hexdigest()
    return hashes[str(file)]
def check(task,name,actual,expected,files=[],status=None,reason=''):
    if status is None:status='PASS' if actual==expected else 'FAIL'
    checks.append(dict(task_id=task,test=name,environment=env,access_scope=scope,relocated=a.relocated,input_files=json.dumps([str(f) for f in files]),input_hashes=json.dumps({str(f):digest(f) for f in files}),actual=json.dumps(actual,ensure_ascii=False,default=str),expected=json.dumps(expected,ensure_ascii=False,default=str),status=status,failure_reason=reason))
try:
    full=cells(root)
    check('TEST4_UNIFIED','unique_metadata',[len(full),int(full.cell_uid.nunique())],[818148,818148] if public else [2474742,2474742],[metadata_root(root)/'metadata/cell_index.parquet'])
    expectations=pd.read_csv(root/'examples/FROZEN_MARKER_EXPECTATIONS.tsv',sep='\t',na_values=['NULL'],keep_default_na=False)
    for source,task in [('Berg_ALI','TEST1_HUMAN'),('RelA_mouse','TEST2_MOUSE'),('Qin_cotton_rat','TEST3_COTTON')]:
        file=root/'expression'/(source+'.h5ad')
        if not file.is_file():
            check(task,'source_H5AD_availability','not included','licensed source H5AD',status='BLOCKED_BY_RELEASE_SCOPE',reason='RelA_mouse is outside confirmed public expression scope.')
            continue
        e=expectations[expectations.source_label.eq(source)]
        measured=marker_counts(root,source,e.cell_uid.unique(),e.native_gene_id.unique())
        j=measured.merge(e[['cell_uid','native_gene_id','measured','count']],on=['cell_uid','native_gene_id'],validate='one_to_one',suffixes=['_actual','_expected'])
        equal=(j.measured_actual==j.measured_expected)&(j.count_actual.eq(j.count_expected)|(j.count_actual.isna()&j.count_expected.isna()))
        check(task,'predefined_marker_counts',[len(j),int((~equal).sum())],[len(e),0],[file,root/'examples/FROZEN_MARKER_EXPECTATIONS.tsv'])
        measured.to_csv(output/(source+'_ACTUAL_MARKER_COUNTS.tsv'),sep='\t',index=False,na_rep='NULL')
        h=ad.read_h5ad(file,backed='r');m=full[full.study_id.eq(source)]
        wanted={'Berg_ALI':19600,'RelA_mouse':84967,'Qin_cotton_rat':269944}[source]
        check(task,'cell_and_gene_dimensions',[h.n_obs,h.n_vars],[wanted,{'Berg_ALI':33694,'RelA_mouse':19059,'Qin_cotton_rat':21667}[source]],[file])
        check(task,'index_source_original_ID',[int(h.obs_names.nunique()),bool(h.obs.study_id.astype(str).eq(source).all()),bool(m.sort_values('expression_row_index').cell_uid.to_numpy().tolist()==h.obs_names.tolist())],[wanted,True,True],[file])
        species={'Berg_ALI':'human','RelA_mouse':'mouse','Qin_cotton_rat':'cotton'}[source]
        check(task,'native_feature_namespace',[str(h.uns['feature_namespace']),bool(h.var.species_namespace.astype(str).eq(species).all()),bool(np.array_equal(h.var_names.astype(str),h.var.native_gene_id.astype(str)))],[species,True,True],[file])
        with np.load(root/h.uns['measurement_mask_reference']) as mask:
            check(task,'actual_capture_gene_mask',[bool(np.array_equal(mask['genes'],h.var.native_gene_id.astype(str))),int(len(set(h.obs.qc_group_id.astype(str))-set(mask['libraries'])))],[True,0],[file,root/h.uns['measurement_mask_reference']])
        with h5py.File(file,'r') as hf:check(task,'raw_integer_count_type',[str(hf['X/data'].dtype),bool('raw' not in hf),len(h.layers)],['int32',True,0],[file])
        if source=='RelA_mouse':
            check(task,'Flex_probe_semantics',bool(h.obs.count_technology.astype(str).eq('probe_counts').all()),True,[file])
            check(task,'animal_proxy_not_inferred_physical',bool(h.obs.unit_kind.astype(str).str.contains('proxy|replicate',case=False).all()),True,[file])
        if source=='Qin_cotton_rat':
            hm=ad.read_h5ad(root/'expression/Berg_ALI.h5ad',backed='r')
            check(task,'separate_human_cotton_gene_space',int(len(set(h.var.native_gene_id.astype(str))&set(hm.var.native_gene_id.astype(str)))),0,[file,root/'expression/Berg_ALI.h5ad']);hm.file.close()
            check(task,'homology_does_not_replace_native_ID',bool((~h.var.homology_is_native_species_truth).all()),True,[file])
            if not public:
                ff=metadata_root(root)/'features/feature_metadata.tsv';feature=pd.read_csv(ff,sep='\t',na_values=['NULL'],keep_default_na=False)
                cotton=feature[feature.species_namespace.eq('cotton')]
                check(task,'supplementary_homology_metadata_only',bool(cotton.homology_is_native_species_truth.astype(str).str.lower().isin(['false','0','0.0']).all()),True,[ff])
        h.file.close()
    c=comparisons(root);ids=identity(root,'CD1c_DC')
    if isinstance(c,dict):
        check('TEST4_UNIFIED','complete_comparison_query',c['status'],'complete407/45 query',status='BLOCKED_BY_RELEASE_SCOPE',reason=c['reason'])
        check('TEST4_UNIFIED','complete_annotation_evidence_query',ids['status'],'complete project identity/evidence',status='BLOCKED_BY_RELEASE_SCOPE',reason=ids['reason'])
    else:
        focus=comparisons(root,True)
        check('TEST4_UNIFIED','comparison_states',[len(c),len(focus),int(focus.analysis_suitability.eq('CONDITIONAL').sum()),int(focus.analysis_suitability.eq('HOLD_FOR_VERIFICATION').sum())],[407,45,42,3],[metadata_root(root)/'comparison/comparison_eligibility.tsv'])
        check('TEST4_UNIFIED','identity_query',int(ids.cells_n.sum()),40283,[metadata_root(root)/'metadata/resource.sqlite'])
        ids.to_csv(output/'CD1c_IDENTITY_QUERY.tsv',sep='\t',index=False,na_rep='NULL')
        an=metadata_root(root)/'annotation/cell_annotation.parquet';aa=pq.read_table(an,columns=['cell_uid','gate2b_version']).to_pandas()
        check('TEST4_UNIFIED','complete_annotation_version',[len(aa),int(aa.cell_uid.nunique()),bool(aa.gate2b_version.eq('annotation_gate2b_v0.1').all())],[2474742,2474742,True],[an])
    u=full[full.study_id.eq('Berg_ALI')].cell_uid.iloc[0];tr=trace(root,u)
    check('TEST4_UNIFIED','cell_to_source_file_trace',[len(tr),bool((root/tr[0]['source_h5ad_file']).is_file()),bool(int(tr[0]['expression_row_index'])>=0)],[1,True,True],[metadata_root(root)/'metadata/cell_index.parquet'])
    full[['study_id','source_id','source_h5ad_file']].drop_duplicates().to_csv(output/'SOURCE_TO_H5AD_QUERY.tsv',sep='\t',index=False)
except Exception as e:
    check('EXECUTION','exception',type(e).__name__+': '+str(e),'successful execution',status='FAIL',reason=str(e))
    import traceback
    (output/'ERROR.txt').write_text(traceback.format_exc())
finally:
    check('TEST5_PORTABILITY' if a.relocated else 'EXECUTION','forbidden_directory_read_attempts',len(denied),0)
    pd.DataFrame(checks).to_csv(output/'H5AD_PORTABILITY_TESTS.tsv',sep='\t',index=False,na_rep='NULL')
    (output/'AUDIT.json').write_text(json.dumps(dict(denied_attempts=denied,opened_paths=sorted(reads)),indent=2)+'\n')
    counts=pd.Series([x['status'] for x in checks]).value_counts().to_dict()
    (output/'SUMMARY.json').write_text(json.dumps(dict(scope=scope,relocated=a.relocated,environment=json.loads(env),checks=counts,elapsed_seconds=time.perf_counter()-begin,second_physical_lab=False,actual_source_files=True),indent=2)+'\n')
    print(json.dumps(dict(scope=scope,relocated=a.relocated,checks=counts)),flush=True)
if any(x['status']=='FAIL' for x in checks):sys.exit(1)

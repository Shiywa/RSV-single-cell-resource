"""Rebuild frozen retained counts from original source assets, without new QC/model fitting.

All locations are command-line arguments. The supplied frozen selection and feature
maps are provenance inputs; successful reconstruction does not grant release rights.
"""
from pathlib import Path
import sys, argparse, json, gzip, tarfile, zipfile, io, re, hashlib, gc, time
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
threadpool_limits(4)

def digest_array(x):
    h=hashlib.sha256();x=np.asarray(x)
    if x.dtype.kind in 'OUS':h.update('\n'.join(map(str,x)).encode())
    else:h.update(np.ascontiguousarray(x).tobytes())
    return h.hexdigest()

def vector(group,key):
    z=group[key]
    if isinstance(z,h5py.Group):
        if 'codes' in z:
            c=z['codes'][:];v=vector(z,'categories');out=np.asarray(v,dtype=object)[np.maximum(c,0)];out[c<0]='';return out
        if 'values' in z:
            out=vector(z,'values');out[z['mask'][:]]='';return out
        raise ValueError('Unsupported H5 vector: '+z.name)
    if '__categories' in group and key in group['__categories']:
        c=z[:];v=vector(group['__categories'],key);out=np.asarray(v,dtype=object)[np.maximum(c,0)];out[c<0]='';return out
    return z.asstr()[:] if z.dtype.kind in 'OSU' else z[:]

def frozen_rows(h,rows):
    """Read complete sparse values for exact requested native rows."""
    rows=np.asarray(rows,dtype=np.int64);order=np.argsort(rows);sorted_rows=rows[order]
    assert len(np.unique(rows))==len(rows)
    splits=np.r_[0,np.flatnonzero(np.diff(sorted_rows)>64)+1,len(rows)]
    blocks=[];ncols=int(h['X'].attrs['shape'][1])
    for a,b in zip(splits[:-1],splits[1:]):
        wanted=sorted_rows[a:b];start=int(wanted[0]);end=int(wanted[-1])+1
        ptr=h['X/indptr'][start:end+1];lo=int(ptr[0]);hi=int(ptr[-1])
        x=sparse.csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],ptr-lo),shape=(end-start,ncols))
        blocks.append(x[wanted-start])
    joined=sparse.vstack(blocks,format='csr');return joined[np.argsort(order)]

class OriginalSources:
    def __init__(self,args):
        self.args=args;self.archives={};self.temp=args.output/'original_members';self.temp.mkdir(exist_ok=True)
    def member(self,path,name):
        if path not in self.archives:self.archives[path]=tarfile.open(path,'r')
        return self.archives[path].extractfile(name)
    def byte_stream(self,path,name):
        return self.member(path,name) if name else path.open('rb')
    def table(self,path,name):
        with self.byte_stream(path,name) as raw:
            with gzip.GzipFile(fileobj=raw) as z:return pd.read_csv(z,sep='\t',header=None,dtype=str,keep_default_na=False)
    def mtx(self,r):
        path=self.args.source_dir/r['source_relpath'];name=r['archive_member']
        if path.suffix=='.zip':
            with zipfile.ZipFile(path) as z:
                with z.open(name) as raw,gzip.GzipFile(fileobj=raw) as f:x=mmread(f).T.tocsr()
                with z.open(name.replace('matrix.mtx.gz','features.tsv.gz')) as raw,gzip.GzipFile(fileobj=raw) as f:features=pd.read_csv(f,sep='\t',header=None,dtype=str,keep_default_na=False)
                with z.open(name.replace('matrix.mtx.gz','barcodes.tsv.gz')) as raw,gzip.GzipFile(fileobj=raw) as f:bars=pd.read_csv(f,sep='\t',header=None,dtype=str,keep_default_na=False)[0].to_numpy()
        else:
            with self.byte_stream(path,name) as raw,gzip.GzipFile(fileobj=raw) as f:x=mmread(f).T.tocsr()
            if r.get('features_override'):
                features=self.table(self.args.source_dir/r['features_override'],'')
            else:features=self.table(path,name.replace('matrix.mtx.gz','features.tsv.gz')) if name else self.table(Path(str(path).replace('matrix.mtx.gz','features.tsv.gz')),'')
            bars=self.table(path,name.replace('matrix.mtx.gz','barcodes.tsv.gz'))[0].to_numpy() if name else self.table(Path(str(path).replace('matrix.mtx.gz','barcodes.tsv.gz')),'')[0].to_numpy()
        ids=features[0].to_numpy()
        if r['study_id']=='Sunshine':
            host=np.asarray([v.startswith('GRCh38_') for v in ids]);x=x[:,host];ids=np.asarray([v.replace('GRCh38_','',1) for v in ids[host]])
        return x,ids,bars,dict(original_cells_n=len(bars),original_features_n=len(features),source_feature_subset='human features only for Sunshine; otherwise original features')
    def h5(self,r):
        path=self.args.source_dir/r['source_relpath'];name=r['archive_member']
        if r['study_id']=='pathogen_specific_TRM':
            part=re.search(r'HIPC-1_(LG|LN)-(CD4|CD8)',path.name);part=part.group(1)+'-'+part.group(2)
            path=self.args.trm_exports/(part+'.h5');name=''
        elif name:
            dest=self.temp/Path(name).name
            if not dest.exists():
                with self.member(path,name) as raw,dest.open('wb') as out:
                    for b in iter(lambda:raw.read(8*1024**2),b''):out.write(b)
            path=dest
        with h5py.File(path,'r') as h:
            g=h['matrix'];shape=tuple(g['shape'][:]);x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=shape).T.tocsr()
            ids=vector(g['features'],'id');bars=vector(g,'cell_uid') if r['study_id']=='pathogen_specific_TRM' else vector(g,'barcodes')
        original_cells=len(bars)
        if r['study_id']=='pathogen_specific_TRM':
            evidence=json.loads((self.args.trm_exports/(part+'_INPUT.json')).read_text());original_cells=evidence['original_cells']
        return x,ids,bars,dict(original_cells_n=original_cells,original_features_n=len(ids),source_feature_subset='original gene/RNA counts; no rounding of normalized values',source_count_file=str(path))
    def h5ad(self,r):
        path=self.args.source_dir/r['source_relpath']
        with h5py.File(path,'r') as h:
            g=h[r['assay_or_layer']];shape=tuple(g.attrs['shape']);encoding=str(g.attrs.get('encoding-type',''))
            cls=sparse.csc_matrix if 'csc' in encoding else sparse.csr_matrix
            x=cls((g['data'][:],g['indices'][:],g['indptr'][:]),shape=shape).tocsr()
            v=h['raw/var'] if r['study_id']=='HLCA' else h['var'];ids=vector(v,'_index') if r['study_id']=='HLCA' else vector(v,'gene_ids')
            symbols=vector(v,'feature_name') if r['study_id']=='HLCA' else vector(v,v.attrs['_index'])
            obs=h['obs'];key=obs.attrs.get('_index','_index');key=key.decode() if isinstance(key,bytes) else key
            bars=vector(obs,key)
        return x,ids,bars,dict(original_cells_n=len(bars),original_features_n=len(ids),source_feature_subset=r['assay_or_layer']+' exact original count layer',feature_symbols=symbols)
    def csv(self,r):
        path=self.args.source_dir/r['source_relpath'];name=r['archive_member']
        if r['study_id']=='CD1c_DC':
            with self.member(path,name) as raw,gzip.GzipFile(fileobj=raw) as z:
                reader=pd.read_csv(z,comment='#',dtype={'Cell_Index':str},chunksize=256)
                chunks=[];bars=[];ids=None
                for c in reader:
                    if ids is None:ids=[col for col in c.columns if col!='Cell_Index' and '|pAbO' not in col]
                    bars.extend(c.Cell_Index.tolist());chunks.append(sparse.csr_matrix(c[ids].to_numpy(dtype=np.int64)))
            x=sparse.vstack(chunks,format='csr');ids=np.asarray(ids);bars=np.asarray(bars)
        else:
            with self.member(path,name) as raw,gzip.GzipFile(fileobj=raw) as z:df=pd.read_csv(z,sep='\t',index_col=0)
            ids=df.index.astype(str).to_numpy();bars=df.columns.astype(str).to_numpy();x=sparse.csr_matrix(df.to_numpy().T)
        return x,ids,bars,dict(original_cells_n=len(bars),original_features_n=len(ids),source_feature_subset='RNA features excluding AbSeq only for Rhapsody; original read-count features for SORT-seq')
    def load(self,r):
        st=r['study_id'];member=r['archive_member'];path=r['source_relpath']
        if st in ['HLCA','Yoshida_airway']:return self.h5ad(r)
        if st in ['CD1c_DC','nasal_submerged_model']:return self.csv(r)
        if st in ['RSV_history_wheeze','iAT2_macrophage','pathogen_specific_TRM']:return self.h5(r)
        return self.mtx(r)
    def close(self):
        for t in self.archives.values():t.close()

def save_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,default=str)+'\n');tmp.replace(path)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['source-dir','selection','feature-dir','frozen-dir','mask-dir','inventory','output']:
        parser.add_argument('--'+name,required=True,type=Path)
    parser.add_argument('--studies',required=True);parser.add_argument('--trm-exports',type=Path)
    parser.add_argument('--export-studies',default='')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    meta=pd.read_parquet(args.selection);inventory=json.loads(args.inventory.read_text());wanted=args.studies.split(',');export_set=set(args.export_studies.split(','));reader=OriginalSources(args)
    for study in wanted:
        if study=='Berg_ALI':continue
        started=time.time();m=meta[meta.study_id==study].sort_values('native_row_index').reset_index(drop=True);sp=m.species_namespace.iloc[0]
        features=pd.read_csv(args.feature_dir/(sp+'_NATIVE_FEATURE_ORDER.tsv'),sep='\t',dtype={'native_gene_id':str},keep_default_na=False);target=pd.Index(features.native_gene_id)
        map_index=pd.read_csv(args.feature_dir/(sp+'_MAP_INDEX.tsv'),sep='\t').set_index('qc_group_id').relative_file.to_dict();mask=np.load(args.mask_dir/(sp+'_gene_availability.npz'));byqc={str(k):v for k,v in zip(mask['libraries'],mask['availability'])};assert np.array_equal(mask['genes'],target.to_numpy())
        groups=[];covered=[];export_chunks=[];export_uids=[];failures=[]
        print(json.dumps(dict(stage='SOURCE_START',study=study,retained_cells=len(m))),flush=True)
        with h5py.File(args.frozen_dir/(sp+'.h5ad'),'r',rdcc_nbytes=512*1024**2) as frozen:
            for r in [x for x in inventory if x['study_id']==study]:
                try:
                    x,ids,bars,original=reader.load(r);x.sum_duplicates();x.sort_indices()
                    assert x.shape==(len(bars),len(ids));assert not pd.Index(bars).has_duplicates
                    assert np.isfinite(x.data).all() and np.all(x.data>=0)
                    if x.data.dtype.kind not in 'iu':assert np.all(x.data==np.floor(x.data))
                    key=m.cell_uid if study=='pathogen_specific_TRM' else m.original_barcode
                    # A GSM is the original capture except pooled iAT2 and reference objects.
                    candidate=m
                    if study not in ['HLCA','Yoshida_airway','Qin_cotton_rat','iAT2_macrophage','pathogen_specific_TRM']:
                        gsm=re.search(r'GSM\d+',r['archive_member'] or r['source_relpath'])
                        if gsm:candidate=m[m.library_id==gsm.group()]
                    if study=='Qin_cotton_rat':candidate=m[m.library_id==Path(r['source_relpath']).stem]
                    candidate_key=candidate.cell_uid if study=='pathogen_specific_TRM' else candidate.original_barcode
                    idx=pd.Index(bars).get_indexer(candidate_key);present=idx>=0;candidate=candidate[present].copy();idx=idx[present]
                    assert len(candidate)>0,'no frozen cells match original object'
                    for qc,b in candidate.groupby('qc_group_id',sort=False):
                        positions=candidate.index.get_indexer(b.index);rows=idx[positions];selected=x[rows].astype(np.int64)
                        mapping=pd.read_csv(args.feature_dir/map_index[qc],sep='\t',dtype=str,keep_default_na=False)
                        id_collision=mapping.groupby('native_gene_id').canonical_gene_id.nunique().gt(1).any()
                        if id_collision:
                            assert 'feature_symbols' in original,'duplicate original IDs require exact original symbols'
                            assert not mapping.groupby(['native_gene_id','native_symbol']).canonical_gene_id.nunique().gt(1).any()
                            keys=pd.MultiIndex.from_arrays([ids,original['feature_symbols']])
                            lookup=mapping.drop_duplicates(['native_gene_id','native_symbol']).set_index(['native_gene_id','native_symbol']).canonical_gene_id
                            canonical=lookup.reindex(keys).fillna('NOT_IN_FROZEN_SOURCE_MAP').to_numpy()
                        else:
                            lookup=mapping.drop_duplicates('native_gene_id').set_index('native_gene_id').canonical_gene_id
                            canonical=lookup.reindex(ids).fillna('NOT_IN_FROZEN_SOURCE_MAP').to_numpy()
                        cols=target.get_indexer(canonical);valid=cols>=0
                        transform=sparse.csr_matrix((np.ones(valid.sum(),np.int64),(np.flatnonzero(valid),cols[valid])),shape=(len(ids),len(target)))
                        rebuilt=(selected@transform).tocsr();rebuilt.sum_duplicates();rebuilt.sort_indices()
                        actual=frozen_rows(frozen,b.native_row_index.to_numpy()).astype(np.int64)
                        delta=(rebuilt-actual).tocsr();delta.eliminate_zeros()
                        measured=np.zeros(len(target),bool);measured[cols[valid]]=True;mask_differences=int(np.count_nonzero(measured!=byqc[qc]))
                        groups.append(dict(study_id=study,object_id=r['object_id'],qc_group_id=qc,retained_cells_n=len(b),original_source_file=r['source_relpath'],archive_member=r['archive_member'],source_cells_n=original['original_cells_n'],source_features_n=len(ids),mapped_source_features_n=int(valid.sum()),unmapped_source_features_n=int((~valid).sum()),target_features_n=len(target),native_counts_difference_n=int(delta.nnz),gene_mask_difference_n=mask_differences,source_cell_join_exact=True,selected_uid_sha256=digest_array(b.cell_uid.to_numpy()),source_feature_order_sha256=digest_array(ids),reconstructed_CSR_data_sha256=digest_array(rebuilt.data.astype(np.int64)),reconstructed_CSR_indices_sha256=digest_array(rebuilt.indices.astype(np.int64)),reconstructed_CSR_indptr_sha256=digest_array(rebuilt.indptr.astype(np.int64)),technical_status='PASS' if delta.nnz==0 and mask_differences==0 else 'FAIL',permission_status='SEPARATE_ASSET_REVIEW_REQUIRED',QC_selection='frozen retained keys, no QC refit'))
                        covered.extend(b.cell_uid.tolist())
                        if study in export_set:export_chunks.append(rebuilt);export_uids.extend(b.cell_uid.tolist())
                    print(json.dumps(dict(stage='OBJECT_DONE',study=study,object=r['object_id'],selected=len(candidate),groups=len(groups),counts_differences=sum(z['native_counts_difference_n'] for z in groups))),flush=True)
                    pd.DataFrame(groups).to_csv(args.output/(study+'_LIBRARY_TESTS.tsv'),sep='\t',index=False)
                    del x;gc.collect()
                except Exception as e:
                    failures.append(dict(object_id=r['object_id'],error=type(e).__name__+': '+str(e)));print(json.dumps(dict(stage='OBJECT_ERROR',study=study,**failures[-1])),flush=True)
        counts=pd.Series(covered).value_counts();missing=set(m.cell_uid)-set(covered);duplicated=int(counts.gt(1).sum());differences=sum(z['native_counts_difference_n'] for z in groups);mask_diff=sum(z['gene_mask_difference_n'] for z in groups)
        passed=not failures and not missing and not duplicated and differences==0 and mask_diff==0
        result=dict(study_id=study,technical_status='PASS' if passed else 'FAIL',retained_cells_n=len(m),verified_unique_cells_n=len(set(covered)),libraries_verified_n=len(set(g['qc_group_id'] for g in groups)),original_objects_tested_n=len([x for x in inventory if x['study_id']==study]),missing_cells_n=len(missing),duplicated_cells_n=duplicated,native_counts_difference_n=differences,gene_mask_difference_n=mask_diff,failures=failures,cell_uid_selection_sha256=digest_array(m.cell_uid.to_numpy()),elapsed_seconds=time.time()-started,annotation_or_QC_refit=False,raw_source_reconstruction=True,technical_result_not_permission_clearance=True,limitations='Uses frozen retained selection and exact historical source-feature transformation, not a new QC run. Original files previously acquired; no fresh download is claimed. Selection/feature-map/annotation delivery rights remain separately evaluated.')
        if passed and study in export_set:
            full=sparse.vstack(export_chunks,format='csr');order=pd.Index(export_uids).get_indexer(m.cell_uid);assert (order>=0).all();full=full[order].tocsr()
            dest=args.output/(study+'_RECONSTRUCTED_RAW.h5')
            with h5py.File(dest,'w') as h:
                g=h.create_group('X');g.attrs['shape']=np.asarray(full.shape,np.int64)
                for name,a in [('data',full.data),('indices',full.indices),('indptr',full.indptr)]:g.create_dataset(name,data=a,compression='gzip',compression_opts=1)
                h.create_dataset('original_native_row_index',data=m.native_row_index.to_numpy());h.create_dataset('native_gene_id',data=np.asarray(target.to_numpy(),dtype='S'));h.create_dataset('cell_uid',data=np.asarray(m.cell_uid,dtype='S'))
            result['reconstructed_count_export']=str(dest);del full
        save_json(args.output/(study+'_RESULT.json'),result);print(json.dumps(dict(stage='SOURCE_COMPLETE',**result)),flush=True)
        del export_chunks,m;gc.collect()
    reader.close()

if __name__=='__main__':main()

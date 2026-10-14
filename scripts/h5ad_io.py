"""Portable, explicit-path, chunked CSR H5AD writer and exact verifier.

Adapts the verified frozen-row reader used in source reconstruction. No model,
normalization, feature renaming, cross-species merge or network access is used.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import time,hashlib,json,resource
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy import sparse

def vector(g,key):
    z=g[key]
    if isinstance(z,h5py.Group):
        if 'codes' in z:
            c=z['codes'][:];v=vector(z,'categories');out=np.asarray(v,dtype=object)[np.maximum(c,0)];out[c<0]=None;return out
        if 'values' in z:
            out=vector(z,'values');out[z['mask'][:]]=None;return out
        raise ValueError('Unsupported H5 vector: '+z.name)
    return z.asstr()[:] if z.dtype.kind in 'OUS' else z[:]

def frozen_rows(h,rows):
    """Same bounded contiguous-run CSR strategy as the verified prior reader."""
    rows=np.asarray(rows,dtype=np.int64);assert len(rows)>0 and np.all(np.diff(rows)>0)
    split=np.r_[0,np.flatnonzero(np.diff(rows)>64)+1,len(rows)];blocks=[]
    for a,b in zip(split[:-1],split[1:]):
        wanted=rows[a:b];start=int(wanted[0]);end=int(wanted[-1])+1
        ptr=h['X/indptr'][start:end+1];lo=int(ptr[0]);hi=int(ptr[-1])
        x=sparse.csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],ptr-lo),shape=(end-start,int(h['X'].attrs['shape'][1])))
        blocks.append(x[wanted-start])
    return sparse.vstack(blocks,format='csr') if len(blocks)>1 else blocks[0]

def batches(rows,batch_cells=2048):
    for a in range(0,len(rows),batch_cells):yield a,min(a+batch_cells,len(rows)),rows[a:a+batch_cells]

def project(x,columns):
    lookup=np.full(x.shape[1],-1,dtype=np.int32);lookup[columns]=np.arange(len(columns),dtype=np.int32)
    mapped=lookup[x.indices];valid=mapped>=0
    # Source-unmeasured union padding can be removed only if every stored value is zero.
    if not valid.all():
        if np.any(x.data[~valid]!=0):raise ValueError('Measured/nonzero feature outside planned source union; no feature was discarded.')
        out=x[:,columns].tocsr()
    else:out=sparse.csr_matrix((x.data,mapped,x.indptr),shape=(x.shape[0],len(columns)))
    return out

def categorical_frame(frame):
    out=frame.copy()
    for k in out:
        if out[k].dtype.kind in 'OUS':out[k]=pd.Categorical(out[k])
    return out

def exact_verify(native,path,rows,columns,maskpath,expected_uids,expected_original,source_label):
    start=time.perf_counter();counts=dict(value_differences=0,index_differences=0,pointer_differences=0,unmeasured_nonzero_values=0,nonzero_outside_source_union=0,entries_compared=0,rows_compared=0,stored_zero_entries=0)
    stream_hash_in=hashlib.sha256();stream_hash_out=hashlib.sha256()
    with h5py.File(native,'r') as n,h5py.File(path,'r') as h,np.load(maskpath) as m:
        nativegenes=vector(n['var'],'canonical_gene_id');genes=vector(h['var'],'native_gene_id')
        assert np.array_equal(genes,nativegenes[columns])
        assert np.array_equal(m['genes'],genes)
        assert np.array_equal(m['native_feature_indices'],columns)
        assert len(genes)==len(np.unique(genes))
        obskey=h['obs'].attrs['_index'];uids=vector(h['obs'],obskey)
        assert np.array_equal(uids,expected_uids)
        assert np.array_equal(vector(h['obs'],'cell_uid'),expected_uids)
        assert np.array_equal(vector(h['obs'],'original_cell_id'),expected_original)
        assert np.array_equal(vector(h['obs'],'native_row_index'),rows)
        qcs=vector(h['obs'],'qc_group_id');ms={x:a for x,a in zip(m['libraries'],m['availability'])}
        assert h['X/data'].dtype==n['X/data'].dtype
        assert tuple(h['X'].attrs['shape'])==(len(rows),len(columns))
        assert len(h['X/indptr'])==len(rows)+1
        ptr=h['X/indptr'][:];assert ptr[0]==0 and ptr[-1]==len(h['X/data'])==len(h['X/indices']) and np.all(np.diff(ptr)>=0)
        native_obs_key=n['obs'].attrs['_index']
        for a,b,r in batches(rows):
            original_keys=n['obs'][native_obs_key].asstr()[r]
            assert np.array_equal(original_keys,expected_original[a:b]),'Original cell ID/native row mismatch'
            assert np.array_equal(np.asarray(['rsvcell:'+str(h['uns/feature_namespace'].asstr()[()])+':'+x for x in original_keys]),expected_uids[a:b]),'Cell UID derived from frozen row differs'
            raw=frozen_rows(n,r);expected=project(raw,columns)
            lo,hi=int(ptr[a]),int(ptr[b]);data=h['X/data'][lo:hi];index=h['X/indices'][lo:hi];p=ptr[a:b+1]-lo
            actual=sparse.csr_matrix((data,index,p),shape=(b-a,len(columns)))
            assert actual.has_sorted_indices and actual.has_canonical_format
            if len(data)!=len(expected.data):raise ValueError('Sparse structure length differs')
            counts['value_differences']+=int(np.count_nonzero(data!=expected.data))
            counts['index_differences']+=int(np.count_nonzero(index!=expected.indices))
            counts['pointer_differences']+=int(np.count_nonzero(p!=expected.indptr))
            rr=np.repeat(np.arange(b-a),np.diff(p));positive=data!=0
            if positive.any():
                for qc in np.unique(qcs[a:b]):
                    selected=(qcs[a:b][rr]==qc)&positive
                    counts['unmeasured_nonzero_values']+=int((~ms[qc][index[selected]]).sum())
            counts['entries_compared']+=len(data);counts['rows_compared']+=b-a;counts['stored_zero_entries']+=int((data==0).sum())
            for v,e in [(data,expected.data),(index,expected.indices),(p.astype(np.int64),expected.indptr.astype(np.int64))]:stream_hash_out.update(np.ascontiguousarray(v).tobytes());stream_hash_in.update(np.ascontiguousarray(e).tobytes())
        assert all(counts[k]==0 for k in ['value_differences','index_differences','pointer_differences','unmeasured_nonzero_values'])
        assert stream_hash_in.hexdigest()==stream_hash_out.hexdigest()
        # Backed AnnData decoding actually checks format interoperability, not h5py alone.
    begin=time.perf_counter();a=ad.read_h5ad(path,backed='r');loadseconds=time.perf_counter()-begin
    assert a.n_obs==len(rows) and a.n_vars==len(columns) and a.obs_names.is_unique and a.var_names.is_unique
    assert a.raw is None and len(a.layers)==0
    assert a.obs['study_id'].astype(str).eq(source_label).all()
    a.file.close()
    return dict(**counts,logical_input_sha256=stream_hash_in.hexdigest(),logical_output_sha256=stream_hash_out.hexdigest(),cell_uid_unique=True,original_ID_exact=True,native_gene_ID_order_exact=True,mask_native_indices_exact=True,count_dtype_preserved=True,complete_value_verification=True,validation_seconds=time.perf_counter()-start,anndata_backed_read_seconds=loadseconds,peak_RAM_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

def write_stream(native,path,obs,var,columns,maskpath,uns,compression='gzip',level=4):
    """Source-level independent retry: write .partial, verify, then atomic rename."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);partial=path.with_suffix('.partial.h5ad')
    if path.exists():raise FileExistsError('Verified target already exists; use its receipt or explicit new target.')
    if partial.exists():raise FileExistsError('Partial file exists; preserve failed artifact before retry.')
    rows=obs.native_row_index.to_numpy(np.int64);assert np.all(np.diff(rows)>0)
    start=time.perf_counter()
    with h5py.File(native,'r') as h:dtype=h['X/data'].dtype;idtype=h['X/indices'].dtype
    shell=ad.AnnData(X=sparse.csr_matrix((len(obs),len(columns)),dtype=dtype),obs=categorical_frame(obs),var=categorical_frame(var),uns=uns)
    shell.write_h5ad(partial,compression=compression,compression_opts=level if compression=='gzip' else None)
    del shell
    with h5py.File(native,'r') as n,h5py.File(partial,'r+') as h:
        del h['X'];g=h.create_group('X');g.attrs['encoding-type']='csr_matrix';g.attrs['encoding-version']='0.1.0';g.attrs['shape']=np.asarray([len(obs),len(columns)],np.int64)
        kw=dict(compression=compression)
        if compression=='gzip':kw['compression_opts']=level
        for key,dt in [('data',dtype),('indices',idtype)]:g.create_dataset(key,shape=(0,),maxshape=(None,),chunks=(1000000,),dtype=dt,**kw)
        g.create_dataset('indptr',shape=(len(obs)+1,),chunks=(min(100000,len(obs)+1),),dtype=np.int64,**kw);g['indptr'][0]=0
        offset=0
        for a,b,r in batches(rows):
            x=project(frozen_rows(n,r),columns);nnz=len(x.data)
            for key,value in [('data',x.data),('indices',x.indices)]:g[key].resize((offset+nnz,));g[key][offset:offset+nnz]=value
            g['indptr'][a+1:b+1]=x.indptr[1:].astype(np.int64)+offset;offset+=nnz
        h.flush()
    writing=time.perf_counter()-start
    verify=exact_verify(native,partial,rows,columns,maskpath,obs.cell_uid.to_numpy(),obs.original_cell_id.to_numpy(),str(obs.study_id.iloc[0]))
    partial.rename(path)
    return dict(writing_seconds=writing,output_bytes=path.stat().st_size,compression=compression or 'none',gzip_level=level if compression=='gzip' else None,**verify)

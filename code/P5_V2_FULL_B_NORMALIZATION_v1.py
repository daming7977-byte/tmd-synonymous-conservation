#!/usr/bin/env python3
"""Validation 3 only: sum-TPM B, frozen +30..+90 outcomes, unchanged models.
Run: python3 -B P5_V2_FULL_B_NORMALIZATION_v1.py
Writes only beside this script. Imports frozen model functions by AST, never
executes the original scripts' top-level writes or biological adjudication.
"""
from pathlib import Path
import sys, os, ast, hashlib, json, platform, datetime, gzip
sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
V1 = OUT.parent/'01_AUTHOR_REPRODUCTION'
sys.path.insert(0, str(V1/'dependencies'))
sys.path.insert(0, str(OUT/'dependencies'))
import numpy as np
import pandas as pd
import openpyxl
import statsmodels
import statsmodels.api as sm
MPT = ROOT/'13_MPT_OUTCOME_REVEAL'
OST = ROOT/'14_OSTA_SPECIFICITY'
BOOK = ROOT/'01_METADATA_SUPP_AUDIT/PUBLICATION_SOURCE_DATA/41594_2025_1691_MOESM3_ESM.xlsx'
PROFILES = [('TMCO1_rep1',84,83),('TMCO1_rep2',82,81),('CCDC47_rep1',80,79),('CCDC47_rep2',78,77),('Nicalin_rep1',76,75),('Nicalin_rep2',74,73),('OST48_rep1',90,89),('OST48_rep2',88,87),('RPN2_rep1',86,85)]
PREFIX = 'P5_V2_B_'
logfile = open(OUT/(PREFIX+'NORMALIZATION_RUN.log'),'w',buffering=1)
def log(x):
 print(x,flush=True); print(x,file=logfile,flush=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
 return h.hexdigest()
inputs={}
def track(p):
 p=Path(p); inputs[p]=sha(p); return p

def read(p): return pd.read_csv(track(p),sep='\t',float_precision='round_trip')
def write(df,name):
 df.to_csv(OUT/name,sep='\t',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)

def main():
 log('START_UTC = '+datetime.datetime.now(datetime.timezone.utc).isoformat())
 log('COMMAND = '+sys.executable+' -B '+str(Path(__file__).resolve()))
 log('SOFTWARE = '+json.dumps(dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,statsmodels=statsmodels.__version__,openpyxl=openpyxl.__version__)))
 track(V1/'SOFTWARE_MANIFEST_SHA256.txt')
 for p in [OUT.parent/'P5_V2_PUBLICATION_VALIDATION_CHARTER.md',V1/'P5_V2_AUTHOR_REPRO_METHOD_ADJUDICATION.md',V1/'P5_V2_AUTHOR_REPRO_SUMMARY.md',ROOT/'07_PREREVEAL_FREEZE/P5_V2_PREREVEAL_MPT_MODEL_SPEC.yaml',OST/'P5_V2_OSTA_SPECIFICITY_FREEZE.yaml']:
  track(p)
 assert inputs[ROOT/'07_PREREVEAL_FREEZE/P5_V2_PREREVEAL_MPT_MODEL_SPEC.yaml']=='181ed0bef6bc81e4750f421c1e0b3f79a7bfd5fc55432d60664a191d7854b912'
 wb=openpyxl.load_workbook(track(BOOK),read_only=True,data_only=True)
 it=wb['DeepTMHMM topology'].iter_rows(values_only=True)
 header=next(it); log('TOPOLOGY_HEADERS = '+repr(header[:5]))
 gi=header.index('Gene name'); ci=header.index('DeepTMHMM protein class')
 classes={}; bg=set(); nrows=0
 for row in it:
  if row[gi]:
   classes.setdefault(row[gi],set()).add(row[ci])
   if row[ci]=='CYTONUC': bg.add(row[gi]); nrows+=1
 conflicts={g:sorted(classes[g]) for g in bg if len(classes[g])>1}
 wb.close()
 log('BACKGROUND_RULE = unique Gene from rows with exact DeepTMHMM protein class == CYTONUC; zeros retained; no paired-positive filter')
 log('CYTONUC_ROWS = '+str(nrows)); log('CYTONUC_BACKGROUND_GENES = '+str(len(bg)))
 log('CLASS_CONFLICTS_INCLUDED_BY_LITERAL_ROW_RULE = '+repr(conflicts))
 counts=read(V1/'CDS_COUNTS_TPM.tsv.gz')
 assert not counts.duplicated(['run','gene']).any()
 constants=[]
 for k,(profile,inp,ip) in enumerate(PROFILES):
  d={'profile':profile,'complex':'MPT' if k<6 else 'OST-A','aggregation':'sum_gene_TPM','background_genes':len(bg),'CYTONUC_source_rows':nrows}
  for label,num in [('IP',ip),('Input',inp)]:
   run=f'SRR336213{num}'; sample=counts[counts.run==run].set_index('gene')
   assert bg <= set(sample.index)
   rpk=sample['count']/(sample.CDS_union_bp/1000)
   assert np.allclose(sample.TPM,rpk/rpk.sum()*1e6,rtol=1e-12,atol=1e-10)
   vals=sample.loc[sorted(bg),'TPM']
   assert np.isfinite(vals).all() and (vals>=0).all()
   d[label+'_run']=run; d['B_'+label]=float(vals.sum()); d[label+'_positive_background_genes']=int((vals>0).sum())
  assert d['B_IP']>0 and d['B_Input']>0
  d['expected_log2_shift']=float(np.log2(d['B_Input']/d['B_IP']))
  constants.append(d)
 write(pd.DataFrame(constants),PREFIX+'CONSTANTS.tsv')
 # B and membership decisions above precede all outcome reads.
 old=read(MPT/'P5_V2_TMD_OUTCOME_TABLE.tsv')
 ost=read(OST/'P5_V2_OSTA_TMD_OUTCOME_TABLE.tsv')
 keys=['human_ensembl_gene','human_ensembl_transcript','p5_tmd_index']
 assert len(old)==5538 and old[keys].equals(ost[keys])
 for c in ['score_z','primary_MPT_TMD_score']:
  np.testing.assert_allclose(old[c],ost[c],rtol=0,atol=1e-12,equal_nan=True)
 log('MPT_TO_OSTA_TABLE_SERIALIZATION_MAX_SCORE_DELTA = '+str(float(np.nanmax(np.abs(old.score_z-ost.score_z)))))
 new=ost.copy(); rows=[]; diagnostics=[]
 # Independently reconstruct every frozen primary-window codon ratio from caches.
 # Frozen arrays and their ratio/log/mean are float32. Preserve the saved TMD
 # statistic exactly and add its B offset in float64; independently evaluate the
 # full codon formula in float64 to quantify reduction/rounding differences.
 for d in constants:
  profile=d['profile']; area=MPT/'occupancy_npz' if d['complex']=='MPT' else OST/'osta_occupancy_npz'
  with np.load(track(area/(d['Input_run']+'.target_codon_rpm.npz')),allow_pickle=False) as z: inp={k:z[k] for k in z.files}
  with np.load(track(area/(d['IP_run']+'.target_codon_rpm.npz')),allow_pickle=False) as z: ip={k:z[k] for k in z.files}
  col=profile+'_TMD_score'; shift=d['expected_log2_shift']
  new[col]=ost[col]+shift
  old_repro=[]; full_direct=[]; codon_res=[]; valid=0
  for idx,r in enumerate(ost.itertuples(index=False)):
   a=inp[r.author_tx]; b=ip[r.author_tx]; A=float(np.mean(a)); end=int(r.author_end_1based)
   saved=float(ost.iloc[idx][col]); direct=np.nan
   if A>0 and end+90<=len(a):
    iw=a[end+29:end+90]; pw=b[end+29:end+90]; assert len(iw)==61
    ratio=(pw+A)/(iw+A); base=np.log2(ratio); rec=float(np.mean(base))
    assert np.isfinite(saved)
    old_repro.append(abs(rec-saved))
    # Promoting after the frozen pseudocount additions preserves their semantics.
    numerator=(pw+A).astype(np.float64); denominator=(iw+A).astype(np.float64)
    full=np.log2((numerator/d['B_IP'])/(denominator/d['B_Input']))
    baseline64=np.log2(numerator/denominator)
    codon_res.append(float(np.max(np.abs(full-baseline64-shift))))
    direct=float(np.mean(full)); full_direct.append(abs(direct-(saved+shift))); valid+=1
   else: assert np.isnan(saved)
   rows.append(dict(record_type='TMD_profile',human_ensembl_gene=r.human_ensembl_gene,human_ensembl_transcript=r.human_ensembl_transcript,p5_tmd_index=r.p5_tmd_index,author_tx=r.author_tx,profile=profile,window='+30..+90',valid=np.isfinite(saved),B_free_log2=saved,full_B_log2=saved+shift,expected_log2_shift=shift,full_B_geometric_E=2**(saved+shift),direct_codon_full_B_mean_log2=direct))
  residual=float(np.nanmax(np.abs(new[col]-ost[col]-shift)))
  assert max(old_repro)<1e-12
  assert max(full_direct)<2e-6 and max(codon_res)<1e-12 and residual<1e-12
  diagnostics.append(dict(kind='profile',analysis=profile,N=valid,expected_constant_shift=shift,max_nonconstant_residual=residual,max_frozen_reconstruction_error=max(old_repro),max_direct_codon_TMD_rounding_difference=max(full_direct),max_codon_algebra_residual=max(codon_res)))
  log('PROFILE = '+json.dumps(diagnostics[-1]))
 for complex_name,col,cs in [('MPT','primary_MPT_TMD_score',constants[:6]),('OST-A','OSTA_TMD_score',constants[6:])]:
  cols=[d['profile']+'_TMD_score' for d in cs]
  new[col]=new[cols].mean(axis=1).where(new[cols].notna().all(axis=1))
  assert new[col].notna().equals(ost[col].notna())
  shift=float(np.mean([d['expected_log2_shift'] for d in cs]))
  residual=float(np.nanmax(np.abs(new[col]-ost[col]-shift)))
  assert residual<1e-12
  diagnostics.append(dict(kind='combined',analysis=complex_name,N=int(new[col].notna().sum()),expected_constant_shift=shift,max_nonconstant_residual=residual))
  for idx,r in enumerate(ost.itertuples(index=False)):
   before=float(ost.iloc[idx][col]); after=float(new.iloc[idx][col])
   rows.append(dict(record_type='TMD_combined',human_ensembl_gene=r.human_ensembl_gene,human_ensembl_transcript=r.human_ensembl_transcript,p5_tmd_index=r.p5_tmd_index,author_tx=r.author_tx,profile=complex_name,window='+30..+90',valid=np.isfinite(before),B_free_log2=before,full_B_log2=after,expected_log2_shift=shift,full_B_geometric_E=2**after))
 write(pd.DataFrame(rows),PREFIX+'NORMALIZED_OUTCOME.tsv.gz')
 # Execute only original frozen function ASTs and literal covariate list.
 code=track(MPT/'P5_V2_BUILD_AND_REVEAL.py').read_text(); tree=ast.parse(code)
 cov=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='continuous_covariates_raw' for t in n.targets))
 env={'np':np,'pd':pd,'sm':sm,'continuous_covariates_z':[c+'_z' for c in cov]}
 nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['cluster_ols','within_gene_model']]
 assert len(nodes)==2
 exec(compile(ast.Module(body=nodes,type_ignores=[]),str(MPT/'P5_V2_BUILD_AND_REVEAL.py'),'exec'),env)
 track(OST/'P5_V2_BUILD_OSTA_SPECIFICITY_v1.py')
 cluster=env['cluster_ols']; within=env['within_gene_model']
 adjusted=['score_z']+env['continuous_covariates_z']+['downstream_segment_lumenal_indicator']
 def fits(df):
  shared=df[df.primary_MPT_TMD_score.notna() & df.OSTA_TMD_score.notna()].copy()
  shared['difference']=shared.primary_MPT_TMD_score-shared.OSTA_TMD_score
  assert len(shared)==2575
  return {'PRIMARY':cluster(df,'primary_MPT_TMD_score',['score_z']),'ADJUSTED':cluster(df,'primary_MPT_TMD_score',adjusted),'WITHIN_PROTEIN':within(df),'MPT_SHARED':cluster(shared,'primary_MPT_TMD_score',['score_z']),'OSTA_SHARED':cluster(shared,'OSTA_TMD_score',['score_z']),'MPT_MINUS_OSTA':cluster(shared,'difference',['score_z'])}
 before=fits(ost); after=fits(new)
 archived=pd.concat([read(MPT/'P5_V2_REVEAL_CORE_RESULTS.tsv'),read(OST/'P5_V2_OSTA_SPECIFICITY_CORE_RESULTS.tsv')]).set_index('analysis')
 for name,b in before.items():
  a=after[name]; row={'kind':'association','analysis':name,'N':b['N'],'clusters':b['clusters']}
  assert b['N']==a['N']==archived.loc[name,'N'] and b['clusters']==a['clusters']
  for stat in ['beta','se','p']:
   row[stat+'_before']=b[stat];row[stat+'_after']=a[stat];row[stat+'_delta']=a[stat]-b[stat]
   row[stat+'_before_vs_archived_delta']=b[stat]-archived.loc[name,stat]
   assert abs(row[stat+'_delta'])<1e-12
   assert abs(row[stat+'_before_vs_archived_delta'])<1e-10
  diagnostics.append(row);log('MODEL = '+json.dumps(row))
 write(pd.DataFrame(diagnostics),PREFIX+'INVARIANCE_RESULTS.tsv')
 maxres=max(d.get('max_nonconstant_residual',0) for d in diagnostics)
 terminal={'CYTONUC_BACKGROUND_GENES':len(bg),'MPT_PROFILES_WITH_B':'6/6','OSTA_PROFILES_WITH_B':'3/3','MAX_NONCONSTANT_RESIDUAL':maxres,'PRIMARY_BETA_DELTA':after['PRIMARY']['beta']-before['PRIMARY']['beta'],'PRIMARY_P_DELTA':after['PRIMARY']['p']-before['PRIMARY']['p'],'ADJUSTED_BETA_DELTA':after['ADJUSTED']['beta']-before['ADJUSTED']['beta'],'WITHIN_BETA_DELTA':after['WITHIN_PROTEIN']['beta']-before['WITHIN_PROTEIN']['beta'],'OSTA_SHARED_BETA_DELTA':after['OSTA_SHARED']['beta']-before['OSTA_SHARED']['beta'],'MPT_MINUS_OSTA_BETA_DELTA':after['MPT_MINUS_OSTA']['beta']-before['MPT_MINUS_OSTA']['beta'],'B_NORMALIZATION_INVARIANCE_CONFIRMED':'TRUE','B_NORMALIZATION_STATUS':'COMPLETE_WITH_DOCUMENTED_METHOD_FIDELITY_LIMITATIONS','FORMAL_PRIMARY_ADJUDICATION_REMAINS':'INTERMEDIATE_BRIDGE_RESULT','NEXT_ACTION':'STOP_AFTER_VALIDATION_3'}
 summary='''# P5-v2 Validation 3 — full B normalization

Validation 3 only. Biological hypothesis testing remains closed. Formal primary adjudication remains **INTERMEDIATE_BRIDGE_RESULT**; OST-A specificity remains NOT SUPPORTED. No population, window, predictor, covariate, model, or adjudication was selected or retuned.

## Method decision

The [publication Methods](https://www.nature.com/articles/s41594-025-01691-6), “Positional enrichment analyses across transcripts,” explicitly call B “sum expression of background cytonuclear genes”. Accordingly B is the sum of Validation-1 uncentered gene TPM in each sample. It is not the median of gene log ratios or median gene expression used in the separate gene-level analysis. A mean over the identical complete gene set would give the same B ratio but different absolute B values; no such alternative was run.

CYTONUC membership is obtained directly from the Source Data topology sheet by exact class equality and deduplicated Gene. Zeros remain in the sums; no paired-positive or detected-gene filter applies. NPIPA2 has both CYTONUC and SP+TM rows and is included by the requested literal row predicate. Validation 1 excluded this conflict for its median centering, but its gene counts/TPMs are unchanged and reused here. This is an explicit distinction, not an alternate-background optimization. All selected genes have Validation-1 TPM entries in all runs.

The source export contains centered normalized expression, not documented raw TPM or original B constants; therefore it cannot independently identify historical B. Positional Methods mix RPM and TPM terminology, and do not completely specify the background expression units/eligibility. Sum of reconstructed TPM is the most directly supported expression interpretation. There is no demonstrated contradiction in the B aggregation itself. Historical counting annotation and zero handling remain Validation-1 fidelity limitations. This validates adding the supported B term to frozen P5 outcomes, not exact identity to historical author positional values.

## Frozen arithmetic and precision

For every profile C=log2(B_Input/B_IP), so log2(E_full)=log2(E_free)+C. The primary TMD statistic remains mean log2 enrichment over exactly 61 codons at +30..+90 after the author TMD end. Every valid profile window was independently read from frozen codon RPM caches; paired-input A is the original float32 full-CDS mean and has not been changed. No alternate window was computed.

Canonical full-B TMD values add C in float64 to the saved B-free statistic, preserving its original float32 reduction. An independent direct evaluation of ((IP+A)/B_IP)/((Input+A)/B_Input) across all valid window codons was averaged in float64 and recorded separately. Small direct-versus-canonical discrepancies are float32 ratio/log/reduction rounding, not biological changes. Both diagnostics are explicitly reported; max_nonconstant_residual refers to canonical full-B minus saved B-free minus C. The output contains all 5,538 TMDs per profile and combined outcome, retains invalid outcomes as missing, and includes direct codon-derived means and geometric enrichment (2**mean log2 E). It is not an arithmetic mean of E.

MPT combines six profiles; OST-A combines three with equal profile weights and all-profile validity. The 2,796 primary-valid TMDs and 2,575 shared TMDs remain fixed. Frozen score/covariate z values were loaded without recalculation. The original cluster-robust OLS and within-protein function definitions were reused by AST extraction, with no original top-level script execution. Before-fit beta/SE/p were checked against archived results. Within-protein eligibility, demeaning, nuisance columns, cluster correction and t inference remain unchanged.

## Results

'''
 summary+='CYTONUC unique genes: '+str(len(bg))+'; exact CYTONUC rows: '+str(nrows)+'.\n\n'
 summary+='| Profile | B IP (TPM sum) | B Input (TPM sum) | Log2 shift |\n|---|---:|---:|---:|\n'
 for d in constants: summary+=f"| {d['profile']} | {d['B_IP']:.12g} | {d['B_Input']:.12g} | {d['expected_log2_shift']:.12g} |\n"
 summary+='\n| Model | N | Beta before → after | SE before → after | p before → after |\n|---|---:|---|---|---|\n'
 for name,b in before.items():
  a=after[name];summary+=f"| {name} | {b['N']} | {b['beta']:.17g} → {a['beta']:.17g} | {b['se']:.17g} → {a['se']:.17g} | {b['p']:.17g} → {a['p']:.17g} |\n"
 summary+='\nMaximum canonical nonconstant residual: '+str(maxres)+'. Maximum independent direct codon/TMD rounding difference: '+str(max(d.get('max_direct_codon_TMD_rounding_difference',0) for d in diagnostics))+'. Association deltas pass the prespecified 1e-12 numerical tolerance; before versus archived passes 1e-10. Direct codon/TMD rounding tolerance is 2e-6 log2. These are numerical checks only.\n\n'
 summary+='## Reproduction and manifest\n\nRun the command recorded in RUN.log with the bundled Python and read-only Validation-1 numpy/scipy dependencies (their existing software manifest is hashed). statsmodels 0.14.6 and patsy 1.0.1 are isolated in dependencies beside the script; no protected input is modified. SHA256 lists every used input, implementation, dependency file and delivered output, excluding itself. All input hashes are rechecked at completion. The TSV results contain per-profile and combined shifts/residuals and all model before/after beta, SE and p. No codon significance or new biological test is performed.\n\n```text\n'+'\n'.join(f'{k} = {v}' for k,v in terminal.items())+'\n```\n'
 (OUT/(PREFIX+'NORMALIZATION_SUMMARY.md')).write_text(summary)
 for p,h in inputs.items(): assert sha(p)==h, f'INPUT_CHANGED: {p}'
 log('INPUT_HASH_RECHECK = PASS')
 for k,v in terminal.items():log(f'{k} = {v}')
 logfile.close()
 outputs=[p for p in OUT.rglob('*') if p.is_file() and p.name!=PREFIX+'NORMALIZATION_SHA256.txt']
 manifest=['# SHA256: absolute read-only inputs, then implementation/dependencies/outputs; manifest excludes itself.']
 manifest += [f'{h}  {p}' for p,h in sorted(inputs.items())]
 manifest += [f'{sha(p)}  {p}' for p in sorted(outputs)]
 (OUT/(PREFIX+'NORMALIZATION_SHA256.txt')).write_text('\n'.join(manifest)+'\n')

if __name__=='__main__': main()

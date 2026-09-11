#!/usr/bin/env python3
"""Validation 2 only: frozen primary eligibility and descriptive attrition.
No enrichment or association values are parsed. No alternate windows evaluated.
Bundled Python + NumPy; every write is beside this script.
"""
from pathlib import Path
import csv, collections, contextlib, datetime, hashlib, io, json, math, sys, zipfile
import numpy as np
OUT=Path(__file__).resolve().parent
ROOT=OUT.parent.parent
TABLE=ROOT/'13_MPT_OUTCOME_REVEAL/P5_V2_TMD_OUTCOME_TABLE.tsv'
CODE=ROOT/'13_MPT_OUTCOME_REVEAL/P5_V2_BUILD_AND_REVEAL.py'
CHARTER=OUT.parent/'P5_V2_PUBLICATION_VALIDATION_CHARTER.md'
OCC=ROOT/'13_MPT_OUTCOME_REVEAL/occupancy_npz'
PAIRS=[('TMCO1_rep1',84,83),('TMCO1_rep2',82,81),('CCDC47_rep1',80,79),('CCDC47_rep2',78,77),('Nicalin_rep1',76,75),('Nicalin_rep2',74,73)]
TECH=['human_ensembl_gene','human_ensembl_transcript','p5_tmd_index','author_tx','author_end_1based','is_frozen_scored_internal','primary_profiles_valid']
FEATURES=[('frozen_evolutionary_score','POST30_90_4D_EXCESS_CONSERVATION_Z_MEDIAN'),('median_window_aa_identity','median_window_aa_identity'),('human_window_GC3','human_window_GC3'),('human_window_CpG_density','human_window_CpG_density'),('TMD_length_aa','p5_length_aa'),('TMD_hydropathy','TMD_mean_Kyte_Doolittle_hydropathy'),('TMD_index','p5_tmd_index'),('downstream_nonTMD_segment_length_aa','downstream_nonTMD_segment_length_aa'),('downstream_lumenal_indicator','downstream_segment_lumenal_indicator')]

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def tsv(name,rows):
 rows=list(rows)
 with open(OUT/name,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def selected_columns(cols):
 # Whole source is hashed, but only these authorized fields are retained/parsed.
 with open(TABLE) as f:
  reader=csv.reader(f,delimiter='\t');h=next(reader);idx=[h.index(c) for c in cols]
  return [dict(zip(cols,(r[i] for i in idx))) for r in reader]
def header_shapes(path):
 # Read NPY headers only. No IP occupancy data payload or enrichment values read.
 result={}
 with zipfile.ZipFile(path) as z:
  for name in z.namelist():
   assert name.endswith('.npy')
   with z.open(name) as f:
    version=np.lib.format.read_magic(f)
    if version==(1,0):shape,fortran,dtype=np.lib.format.read_array_header_1_0(f)
    elif version==(2,0):shape,fortran,dtype=np.lib.format.read_array_header_2_0(f)
    else:raise ValueError(('Unsupported NPY header',version))
   assert len(shape)==1 and not dtype.hasobject
   result[name[:-4]]=shape[0]
 return result

def main():
 print('START_UTC',datetime.datetime.now(datetime.timezone.utc).isoformat())
 print('COMMAND',sys.executable,str(Path(__file__).resolve()))
 print('SOFTWARE',sys.version.replace('\n',' '),'NumPy',np.__version__)
 print('SCOPE Validation 2 only; technical predicates fixed before feature parsing')
 inputs=[TABLE,CODE,CHARTER]+[OCC/f'SRR336213{n}.target_codon_rpm.npz' for n in range(73,85)]
 hashes={str(p):sha(p) for p in inputs}
 rows=selected_columns(TECH);N=len(rows);genes=np.array([r['human_ensembl_gene'] for r in rows]);tx=[r['author_tx'] for r in rows]
 assert N==5538 and len(set(genes))==1400
 assert all(r['is_frozen_scored_internal']=='1' for r in rows)
 assert len({(r['human_ensembl_gene'],r['p5_tmd_index']) for r in rows})==N
 gene_tx=collections.defaultdict(set)
 for r in rows:gene_tx[r['human_ensembl_gene']].add(r['author_tx'])
 assert all(len(x)==1 for x in gene_tx.values()) and len(set(tx))==1400
 shapes={};baselines={};profile_geometry=[];profile_positive=[];runrows=[]
 for name,inp,ip in PAIRS:
  ipshapes=header_shapes(OCC/f'SRR336213{ip}.target_codon_rpm.npz');inshapes={};means={}
  with np.load(OCC/f'SRR336213{inp}.target_codon_rpm.npz',allow_pickle=False) as z:
   assert set(tx)<=set(z.files)
   for t in sorted(set(tx)):
    a=z[t];assert a.ndim==1 and len(a)>0 and np.all(np.isfinite(a)) and np.all(a>=0)
    inshapes[t]=len(a);means[t]=float(np.mean(a)) # same float32 mean as frozen implementation
    assert (means[t]>0)==bool(np.any(a>0))
  assert all(inshapes[t]==ipshapes[t] for t in set(tx))
  shapes[inp]=inshapes;shapes[ip]=ipshapes;baselines[inp]=means
  geometry=[]
  for r in rows:
   end=int(r['author_end_1based']);s=end+29;e=end+90;t=r['author_tx']
   geometry.append(s>=0 and e-s==61 and e<=inshapes[t] and e<=ipshapes[t])
  geom=np.array(geometry,bool);pos=np.array([means[t]>0 for t in tx]);profile_geometry.append(geom);profile_positive.append(pos)
  runrows.append({'profile':name,'input_run':f'SRR336213{inp}','IP_run_header_only':f'SRR336213{ip}','geometry_valid_TMD_N':int(geom.sum()),'input_A_positive_TMD_N':int(pos.sum()),'input_A_positive_gene_N':len(set(genes[pos]))})
 geometry=np.all(profile_geometry,axis=0);positive=np.all(profile_positive,axis=0);included=geometry&positive
 observed=np.array([int(r['primary_profiles_valid']) for r in rows]);reconstructed=np.sum(np.array(profile_geometry)&np.array(profile_positive),axis=0)
 assert np.array_equal(observed,reconstructed),'Per-row valid-profile count differs from input-only reconstruction'
 assert geometry.sum()==5538 and positive.sum()==2796 and included.sum()==2796
 assert np.array_equal(included,observed==6)
 included.setflags(write=False)
 print('ELIGIBILITY VERIFIED before attributes read: every row matches frozen valid-profile count; no enrichment magnitude or evolutionary value used')
 flows=[];priorN=N;priorG=len(set(genes))
 for order,(stage,label,mask) in enumerate([('prelocked','Prelocked scored internal TMDs',np.ones(N,bool)),('primary_geometry','Complete primary +30 to +90 amino-acid window',geometry),('all_six_input_A_positive','Positive transcript-wide input A_j in all six MPT profiles',geometry&positive),('formal_primary_valid','Frozen formal primary-valid set',observed==6)],1):
  n=int(mask.sum());g=len(set(genes[mask]));flows.append({'stage_order':order,'stage':stage,'manuscript_stage':label,'TMD_N':n,'gene_N':g,'TMD_excluded_from_prior':priorN-n,'gene_excluded_from_prior':priorG-g,'TMD_fraction_retained_from_prior':n/priorN,'TMD_fraction_retained_from_initial_5538':n/N,'gene_fraction_retained_from_prior':g/priorG,'gene_fraction_retained_from_initial_1400':g/1400});priorN=n;priorG=g
 tsv('P5_V2_ATTRITION_FLOW.tsv',flows)
 tsv('PRIMARY_PROFILE_INPUT_A_AUDIT.tsv',runrows)
 tmdrows=[]
 for j,r in enumerate(rows):
  entry={k:r[k] for k in TECH};entry.update(primary_geometry_valid=int(geometry[j]),all_six_input_A_positive=int(positive[j]),reconstructed_profiles_valid=int(reconstructed[j]),formal_primary_valid=int(included[j]))
  for name,inp,ip in PAIRS:entry[name+'_input_A_j']=baselines[inp][r['author_tx']]
  tmdrows.append(entry)
 tsv('PRIMARY_ELIGIBILITY_AUDIT.tsv',tmdrows)
 # Outcomes are never parsed. Frozen score is read now, only as requested diagnostic.
 feature_rows=selected_columns(['human_ensembl_gene','p5_tmd_index']+[c for _,c in FEATURES if c!='p5_tmd_index'])
 assert [(r['human_ensembl_gene'],r['p5_tmd_index']) for r in feature_rows]==[(r['human_ensembl_gene'],r['p5_tmd_index']) for r in rows]
 descriptives=[]
 for label,col in FEATURES:
  def numeric(v):
   try:return float(v)
   except (ValueError,TypeError):return float('nan')
  vals=np.array([numeric(r[col]) for r in feature_rows]);good=np.isfinite(vals);a=vals[included&good];b=vals[~included&good]
  d={'attribute':label,'source_column':col,'included_total_N':int(included.sum()),'excluded_total_N':int((~included).sum())}
  for group,x,mask in [('included',a,included),('excluded',b,~included)]:
   q=np.quantile(x,[.25,.5,.75],method='linear');d.update({group+'_finite_N':len(x),group+'_missing_N':int((mask&~good).sum()),group+'_mean':float(np.mean(x)),group+'_SD_sample':float(np.std(x,ddof=1)),group+'_Q1':float(q[0]),group+'_median':float(q[1]),group+'_Q3':float(q[2]),group+'_IQR_width':float(q[2]-q[0])})
  if label=='downstream_lumenal_indicator':
   assert set(vals[good])<={0.,1.};pa=float(a.mean());pb=float(b.mean());den=math.sqrt((pa*(1-pa)+pb*(1-pb))/2);rule='binary: sqrt((p_I*(1-p_I)+p_E*(1-p_E))/2)'
  else:den=math.sqrt((np.var(a,ddof=1)+np.var(b,ddof=1))/2);rule='continuous/ordinal: sqrt((sample_var_I+sample_var_E)/2)'
  d['SMD_included_minus_excluded']=float((a.mean()-b.mean())/den) if den>0 else float('nan');d['SMD_denominator_rule']=rule;descriptives.append(d)
 tsv('P5_V2_INCLUDED_EXCLUDED_DESCRIPTIVES.tsv',descriptives)
 summary=[]
 for g in sorted(set(genes)):
  ix=np.where(genes==g)[0];n=len(ix);v=int(included[ix].sum());summary.append({'human_ensembl_gene':g,'author_tx':next(iter(gene_tx[g])),'prelocked_internal_TMD_N':n,'primary_geometry_valid_TMD_N':int(geometry[ix].sum()),'formal_valid_TMD_N':v,'excluded_TMD_N':n-v,'fraction_TMD_retained':v/n,'gene_status':'fully_retained' if v==n else 'fully_excluded' if v==0 else 'partially_retained','all_six_input_A_positive':int(positive[ix[0]])})
 tsv('P5_V2_ATTRITION_GENE_SUMMARY.tsv',summary)
 alln=np.array([g['prelocked_internal_TMD_N'] for g in summary]);validn=np.array([g['formal_valid_TMD_N'] for g in summary]);lost=alln-validn
 contributing=validn>0
 assert contributing.sum()==736 and np.sum((validn>0)&(lost>0))==0
 distributions=[];dists=[('all_prelocked_genes_prelocked_TMDs_per_gene',alln,np.ones(1400,int),'gene'),('contributing_genes_prelocked_TMDs_per_gene',alln[contributing],np.ones(int(contributing.sum()),int),'gene'),('formal_valid_TMDs_per_contributing_gene',validn[contributing],np.ones(int(contributing.sum()),int),'gene'),('formal_valid_TMD_weighted_gene_multiplicity',validn[contributing],validn[contributing],'formal_valid_TMD')]
 diststats=[]
 for name,vals,weights,unit in dists:
  for k in sorted(set(vals)):
   w=int(weights[vals==k].sum());distributions.append({'distribution':name,'internal_TMDs_per_gene':int(k),'frequency':w,'unit':unit,'denominator':int(weights.sum()),'fraction':w/weights.sum()})
  expanded=np.repeat(vals,weights);q=np.quantile(expanded,[.25,.5,.75]);diststats.append({'distribution':name,'N':len(expanded),'unit':unit,'min':int(expanded.min()),'Q1':float(q[0]),'median':float(q[1]),'Q3':float(q[2]),'max':int(expanded.max()),'mean':float(expanded.mean())})
 tsv('INTERNAL_TMD_MULTIPLICITY_DISTRIBUTIONS.tsv',distributions)
 tsv('INTERNAL_TMD_MULTIPLICITY_SUMMARY.tsv',diststats)
 loss_sort=np.sort(lost)[::-1];total=int(lost.sum());affected=int((lost>0).sum());cum=np.cumsum(loss_sort)
 half=int(np.searchsorted(cum,.5*total)+1);eighty=int(np.searchsorted(cum,.8*total)+1)
 concentration={'excluded_TMD_N':total,'affected_genes':affected,'affected_gene_fraction':affected/1400,'partially_retained_genes':int(((validn>0)&(lost>0)).sum()),'largest_gene_excluded_TMD_N':int(loss_sort[0]),'largest_gene_share_of_excluded_TMDs':float(loss_sort[0]/total),'genes_needed_for_50pct_losses':half,'genes_needed_for_80pct_losses':eighty,'effective_loss_contributing_genes_inverse_HHI':float(1/np.sum((lost/total)**2))}
 for frac in [.01,.05,.1]:
  k=math.ceil(1400*frac);concentration[f'top_{int(100*frac)}pct_all_genes_N']=k;concentration[f'top_{int(100*frac)}pct_all_genes_loss_share']=float(loss_sort[:k].sum()/total)
 (OUT/'ATTRITION_CONCENTRATION.json').write_text(json.dumps(concentration,indent=2)+'\n')
 maxd=max(abs(d['SMD_included_minus_excluded']) for d in descriptives if np.isfinite(d['SMD_included_minus_excluded']));maxfeature=max(descriptives,key=lambda d:abs(d['SMD_included_minus_excluded']))['attribute']
 terminal=f'''PRELOCKED_TMD = {N}
PRELOCKED_GENES = 1400
PRIMARY_GEOMETRY_VALID_TMD = {int(geometry.sum())}
ALL_SIX_INPUT_A_POSITIVE_TMD = {int(positive.sum())}
FORMAL_PRIMARY_VALID_TMD = {int(included.sum())}
FORMAL_PRIMARY_VALID_GENES = {int(contributing.sum())}
ATTRITION_CAUSE = At_least_one_of_six_input_transcript_A_j_equals_zero
OUTCOME_DEPENDENT_SELECTION_USED = NO
MAX_ABS_STANDARDIZED_MEAN_DIFFERENCE = {maxd:.6f} ({maxfeature})
ATTRITION_CONCENTRATION = BROADLY_DISTRIBUTED_across_{affected}_genes; {half}_genes_account_for_50pct_of_losses
ATTRITION_AUDIT_STATUS = PASS_VALIDATION_2_COMPLETE
NEXT_ACTION = Stop_for_independent_Validation_2_review
'''
 (OUT/'TERMINAL_SUMMARY.txt').write_text(terminal)
 text=['# P5-v2 Validation 2: attrition audit','',
 'The frozen set contains 5,538 internal TMDs from 1,400 genes. All satisfy the primary +30..+90 geometry. Requiring positive transcript-wide input A_j in every one of the six paired MPT profiles retains 2,796 TMDs (50.4875%) from 736 genes (52.5714%). The other 2,742 TMDs belong to 664 fully excluded genes. There are no partially retained genes.','',
 '## Manuscript-ready flow table','',
 '| Stage | TMDs | Genes | TMD retention from prior | TMD retention from initial |',
 '|---|---:|---:|---:|---:|']
 for r in flows:text.append(f"| {r['manuscript_stage']} | {r['TMD_N']:,} | {r['gene_N']:,} | {100*r['TMD_fraction_retained_from_prior']:.2f}% | {100*r['TMD_fraction_retained_from_initial_5538']:.2f}% |")
 text+=['','The formal-valid row confirms the previous stage; it is not an additional filter. P5_V2_ATTRITION_FLOW.tsv also supplies gene retention fractions and losses at each stage.','',
 '## Exact reconstruction and evidence limits','',
 'The outcome table has primary_profiles_valid but no per-profile A_j or explicit geometry flags. The table alone can establish the initial and final populations, but cannot independently establish their cause. This audit therefore uses retained transcript-codon occupancy caches as supporting technical evidence: six input payloads for transcript-wide A_j and all paired input/IP array lengths for geometry. IP occupancy values, enrichment columns, factor values, negative-control values and association results are not parsed. IP array lengths are read from NPY headers only. Entire files are hashed for integrity without interpreting excluded fields.','',
 'The source implementation (P5_V2_BUILD_AND_REVEAL.py, lines 1293–1365 and 1418–1432) defines A_j as mean input codon RPM, uses author_end_1based with the 61 codons end+30 through end+90 inclusive (Python end+29:end+90), and requires six valid profiles. Each paired input/IP array has equal length. Every primary interval is fully contained. All input arrays are finite and nonnegative. Their original float32 means reproduce A_j>0 exactly, also agreeing with the independent any-positive-codon check. No occupancy, window or threshold is changed.','',
 'For all 5,538 rows, the reconstructed number of geometry-valid, positive-A_j profiles equals the frozen primary_profiles_valid count, including counts 0 through 5 among exclusions. Formal inclusion exactly equals the all-six input predicate. Geometry excludes zero rows. Thus this audit finds no inclusion criterion based on evolutionary-score value or enrichment magnitude. This is verified eligibility logic for the frozen data, not a claim that coverage is biologically random or that retained and excluded proteins have identical properties.','',
 'Eligibility is constructed and locked before reading diagnostic attributes. The only evolutionary value subsequently used is the frozen score requested for descriptive comparison; it never enters a selection predicate. No inferential p-values are computed or reported, no model is fitted, and no alternate window is evaluated.','',
 '## Included versus excluded descriptive attributes','',
 'The unit is the prelocked TMD, with included-minus-excluded SMD. Q1/median/Q3 use linear empirical quantiles; IQR width and missingness are in the TSV. Continuous and ordinal SMDs divide the mean difference by sqrt((sample variance included + sample variance excluded)/2). The binary indicator instead uses sqrt((p_I(1-p_I)+p_E(1-p_E))/2); means are its proportions. These are TMD-weighted descriptions, not independent-gene inference. The original raw attributes are used: frozen score (not score_z), P5 TMD length/index, and untransformed downstream segment length.','',
 '| Attribute | Included median [Q1, Q3] | Excluded median [Q1, Q3] | SMD |',
 '|---|---:|---:|---:|']
 for d in descriptives:text.append(f"| {d['attribute']} | {d['included_median']:.5g} [{d['included_Q1']:.5g}, {d['included_Q3']:.5g}] | {d['excluded_median']:.5g} [{d['excluded_Q1']:.5g}, {d['excluded_Q3']:.5g}] | {d['SMD_included_minus_excluded']:.6f} |")
 text+=[f'','Largest absolute descriptive SMD: '+f'{maxd:.6f} for {maxfeature}. No significance or causal interpretation is attached.','',
 '## Internal TMD multiplicity','',
 'Internal TMD counts refer exclusively to the 5,538 prelocked scored internal TMDs, not all TMDs in each protein. Contributing genes are counted once in gene-level distributions. For completeness, formal-valid multiplicity is reported both per contributing gene and weighted by the number of formal-valid TMDs: the latter describes the multiplicity encountered by a randomly selected formal-valid TMD. Full exact histograms are in INTERNAL_TMD_MULTIPLICITY_DISTRIBUTIONS.tsv.','',
 '| Distribution | Unit | N | Median [Q1, Q3] | Range |',
 '|---|---|---:|---:|---:|']
 for d in diststats:text.append(f"| {d['distribution']} | {d['unit']} | {d['N']} | {d['median']:g} [{d['Q1']:g}, {d['Q3']:g}] | {d['min']}–{d['max']} |")
 text+=['','## Concentration of attrition','',
 f"Losses are distributed across {affected}/1,400 genes ({affected/14:.2f}%), rather than being dominated by a handful of genes. The largest single gene contributes {loss_sort[0]} excluded TMDs ({100*loss_sort[0]/total:.2f}% of losses). Reaching half of all excluded TMDs requires {half} genes; reaching 80% requires {eighty} genes. The top 1%, 5% and 10% of all prelocked genes ranked by loss account for {100*concentration['top_1pct_all_genes_loss_share']:.2f}%, {100*concentration['top_5pct_all_genes_loss_share']:.2f}% and {100*concentration['top_10pct_all_genes_loss_share']:.2f}% of losses, respectively. This is descriptive concentration of a gene-level coverage filter; loss is not dispersed independently among TMDs within a gene. No arbitrary inferential concentration test is used.",
 '', '## Manuscript wording','',
 f'“Of 5,538 prelocked internal transmembrane domains from 1,400 genes, all had complete primary downstream windows (+30 to +90 amino acids). The prespecified requirement for positive transcript-wide input coverage in all six MPT profiles retained 2,796 domains from 736 genes. The remaining 2,742 domains belonged to 664 genes with zero input baseline coverage in at least one profile. No domain was excluded for primary-window geometry or on the basis of its evolutionary score or enrichment value. Included and excluded domains were compared descriptively using outcome-independent sequence and topology attributes.”',
 '', '## Reproducibility','',
 'All writes are confined to 02_ATTRITION_AUDIT. The SHA256 manifest covers every input, the implementation, and outputs except the manifest itself. Input hashes are checked before and after execution. The script reads only named technical and descriptive columns; no protected file is modified. No Validation 1 rerun or Validation 3 work was performed. Stop for independent Validation 2 review.',
 '', '```text',terminal.rstrip(),'```']
 (OUT/'P5_V2_ATTRITION_AUDIT_SUMMARY.md').write_text('\n'.join(text)+'\n')
 for p,h in hashes.items():assert sha(p)==h,'Input changed: '+p
 print('INPUT_SHA256_UNCHANGED = TRUE')
 print('CONCENTRATION',json.dumps(concentration));print('VALID_PROFILE_COUNT_DISTRIBUTION',dict(collections.Counter(observed.tolist())))
 print(terminal,end='');print('END_UTC',datetime.datetime.now(datetime.timezone.utc).isoformat())
 return hashes,terminal

if __name__=='__main__':
 with open(OUT/'P5_V2_ATTRITION_AUDIT_RUN.log','w') as logfile,contextlib.redirect_stdout(logfile):hashes,terminal=main()
 files=[p for p in OUT.iterdir() if p.is_file() and p.name!='P5_V2_ATTRITION_AUDIT_SHA256.txt']
 manifest=''.join(h+'  '+p+'\n' for p,h in sorted(hashes.items()))+''.join(sha(p)+'  '+str(p)+'\n' for p in sorted(files))
 (OUT/'P5_V2_ATTRITION_AUDIT_SHA256.txt').write_text(manifest)
 print(terminal,end='')

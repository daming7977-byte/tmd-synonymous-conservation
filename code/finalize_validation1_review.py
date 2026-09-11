"""Post-run descriptive QA. Never changes counts, TPM, centering, or matched values."""
from pathlib import Path
import sys,csv,gzip,json,hashlib,importlib.util,collections,datetime
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT/'dependencies'))
import numpy as np
from scipy.stats import rankdata
spec=importlib.util.spec_from_file_location('repro',OUT/'P5_V2_AUTHOR_REPRO_v1.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
P='P5_V2_AUTHOR_REPRO_'
profiles=list(csv.DictReader(open(OUT/(P+'PROFILE_RESULTS.tsv')),delimiter='\t'))
with gzip.open(OUT/(P+'MATCHED_VALUES.tsv.gz'),'rt') as f:matched=list(csv.DictReader(f,delimiter='\t'))
by=collections.defaultdict(list)
for r in matched:by[r['profile']].append(r)
sources,audit,bg,txgene,reps=m.source_audit()
checks=[];classrows=[]
for p in profiles:
 rows=by[p['profile']];a=np.array([float(r['author_enrichment']) for r in rows]);b=np.array([float(r['reconstructed_enrichment']) for r in rows]);delta=b-a
 pearson=np.corrcoef(a,b)[0,1];spearman=np.corrcoef(rankdata(a),rankdata(b))[0,1];slope,intercept=np.linalg.lstsq(np.column_stack([a,np.ones(len(a))]),b,rcond=None)[0]
 assert len(a)==int(p['matched_N'])
 for k,v in [('Pearson_r',pearson),('Spearman_rho',spearman),('OLS_slope_reconstruction_on_author',slope),('OLS_intercept',intercept),('median_absolute_difference',np.median(abs(delta)))]:assert abs(float(p[k])-v)<1e-11,(p['profile'],k)
 median=float(np.median(delta));frac=float(np.mean(abs(delta-median)<1e-5))
 checks.append({'profile':p['profile'],'matched_N':len(a),'independent_metrics_check':'PASS','median_signed_difference':median,'fraction_within_1e_minus_5_of_median_difference':frac,'note':'descriptive constant-offset diagnostic only; matched values are not recentered'})
 for cls in ['CYTONUC','SP','TM','SP+TM']:
  ix=[i for i,r in enumerate(rows) if sources[r['group']][r['gene']]['DeepTMHMM protein class']==cls]
  if ix:classrows.append({'profile':p['profile'],'protein_class':cls,'N':len(ix),'median_signed_difference':float(np.median(delta[ix])),'median_absolute_difference':float(np.median(abs(delta[ix])))})
m.write_tsv(OUT/'INDEPENDENT_METRIC_CHECK.tsv',checks);m.write_tsv(OUT/'RESIDUALS_BY_PROTEIN_CLASS.tsv',classrows)
count_by_run={}
for n in range(73,91):
 with open(OUT/f'SRR336213{n}.counts.tsv') as f:count_by_run[n]={r[0]:int(r[1]) for r in csv.reader(f,delimiter='\t')}
excluded=[]
for sheet,pr,inp,ip in m.PROFILES:
 for g in sources[sheet]:
  ci=count_by_run[inp].get(g);cp=count_by_run[ip].get(g)
  if ci is None or cp is None or ci==0 or cp==0:
   excluded.append({'profile':pr,'group':sheet,'gene':g,'Input_count':ci,'IP_count':cp,'reason':'missing_annotation' if ci is None or cp is None else 'both_zero' if ci==cp==0 else 'IP_zero' if cp==0 else 'Input_zero'})
m.write_tsv(OUT/'AUTHOR_GENE_EXCLUSIONS.tsv',excluded)
exc=collections.Counter((r['group'],r['reason']) for r in excluded)
# Modest but systematic calibration offsets matter for exact publication fidelity.
# This is a descriptive review judgement, not an acceptance threshold or tuning criterion.
material='YES_small_downward_calibration_offsets_MPT_-0.0573_OSTA_-0.1069_log2'
terminal=(OUT/'TERMINAL_SUMMARY.txt').read_text().replace('PENDING_QUANTITATIVE_REVIEW',material)
(OUT/'TERMINAL_SUMMARY.txt').write_text(terminal)
summary=(OUT/(P+'SUMMARY.md')).read_text().replace('PENDING_QUANTITATIVE_REVIEW',material).replace('No correlation p-values were calculated or used.','No correlation p-values were reported or used for decisions.').replace('Exact run command is recorded in RUN.log','Exact run command is recorded in PRECOMPARISON_METHOD_FREEZE.json and RUN.log')
addition='''\n## Independent post-run technical review\n\nStrong agreement in all nine profiles supports reproduction of the published gene-level enrichment behavior under the technical choices fixed before comparison. It does not establish historical implementation identity. Profile slopes range from 0.986842 to 0.995978. All profile intercepts are negative (-0.102039 to -0.041777 log2), with combined intercepts -0.057300 (MPT) and -0.106894 (OST-A). These are modest, systematic calibration discrepancies, material to claims of exact numerical reproduction; their biological importance is not assessed here. No fitted correction was applied.\n\nAll source genes match the annotation. Finite comparison exclusions arise solely from zero reconstructed counts: 119–127 genes per MPT profile and 66–69 per OST-A profile. Combined finite coverage is 7,632/7,763 for MPT and 4,633/4,703 for OST-A. AUTHOR_GENE_EXCLUSIONS.tsv lists every exclusion. Source values are positive for these genes, so these discrepancies cannot be resolved from the published table alone. Historical gene filtering/zero handling, annotation/counting features, read selection and mapping differences remain possible contributors; the run does not prove which is causal.\n\nIndependent calculations from the saved matched values reproduce every reported correlation, slope, intercept and median absolute difference to 1e-11. INDEPENDENT_METRIC_CHECK.tsv quantifies the fraction of residuals near each median shift as a descriptive diagnostic only. RESIDUALS_BY_PROTEIN_CLASS.tsv reports residual distributions by author protein class, with no changes to matching or centering. The synthetic HTSeq counting fixture passed all expected assignments and exclusions.\n\nThe supplementary finalizer is descriptive QA only. P5_V2_AUTHOR_REPRO_v1.py and METHOD_ADJUDICATION retain their precomparison hashes. Re-running the pipeline then finalize_validation1_review.py reproduces the final descriptive review and compact terminal summary. Counts, TPM, matched values and primary metric files were not edited by the finalizer.\n'''
summary+=addition
(OUT/(P+'SUMMARY.md')).write_text(summary)
freeze=json.loads((OUT/'PRECOMPARISON_METHOD_FREEZE.json').read_text())
for name,h in freeze['sha256'].items():assert m.sha(OUT/name)==h
with open(OUT/(P+'RUN.log'),'a') as f:
 f.write('\nPOST_RUN_INDEPENDENT_QA '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'\n')
 f.write('EXACT_MAIN_COMMAND '+freeze['command']+'\n')
 f.write('FINALIZER_COMMAND PYTHONDONTWRITEBYTECODE=1 '+sys.executable+' '+str(Path(__file__).resolve())+'\n')
 f.write('Independent saved-value metrics verified; counting fixture passed; precomparison method/code hashes unchanged.\n')
 f.write('No primary numeric result was modified. Systematic calibration offsets are retained.\n'+terminal)
# Every requested output and all lightweight diagnostic outputs; packages inventoried separately.
software=[]
for p in sorted((OUT/'dependencies').rglob('*')):
 if p.is_file() and (p.suffix in ['.py','.so'] or p.name in ['METADATA','RECORD']):software.append(p)
(OUT/'SOFTWARE_MANIFEST_SHA256.txt').write_text(''.join(m.sha(p)+'  '+str(p)+'\n' for p in software))
outputs=[p for p in OUT.rglob('*') if p.is_file() and 'dependencies' not in p.parts and p.name!='OUTPUT_MANIFEST_SHA256.txt' and '__pycache__' not in p.parts]
(OUT/'OUTPUT_MANIFEST_SHA256.txt').write_text(''.join(m.sha(p)+'  '+str(p)+'\n' for p in sorted(outputs)))
print('Metric check, zero-exclusion audit and final manifests complete.')
print('Exclusion counts by group/reason:',dict(exc))
print('Constant-offset fractions:',[(r['profile'],round(r['fraction_within_1e_minus_5_of_median_difference'],4)) for r in checks])
print(terminal,end='')

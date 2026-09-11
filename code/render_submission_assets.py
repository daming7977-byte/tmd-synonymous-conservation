from pathlib import Path
import csv, gzip, json, re, shutil, hashlib, os, sys
os.environ['MPLCONFIGDIR']=os.environ['P5_MPLCONFIGDIR']
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image, ImageOps, ImageDraw

ROOT=Path(os.environ['P5_PROJECT_ROOT']).resolve()
EVO=ROOT.parent/'04_STAGE1_EVOLUTION'
OUT=Path(os.environ['P5_RENDER_OUT']).resolve()
assert OUT != ROOT and ROOT not in OUT.parents, 'Render outside the frozen project' 
MAN=ROOT/'17_MANUSCRIPT'
OUT.mkdir(exist_ok=True)
for name in ['manuscript','figures','supplement','source_data','reproducibility/code','reproducibility/manifests','reproducibility/results']:
    (OUT/name).mkdir(parents=True,exist_ok=True)
protected={p:p.read_bytes() for p in MAN.glob('*.md')}
sources={}
def source(path):
    p=Path(path)
    if not p.is_absolute():p=ROOT/p
    sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    return p
def rows(path):
    p=source(path)
    op=gzip.open if p.suffix=='.gz' else open
    with op(p,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,text): (OUT/name).write_text(text)
def tsv(name,data):
    keys=list(dict.fromkeys(k for r in data for k in r))
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t');w.writeheader();w.writerows(data)
core=rows('13_MPT_OUTCOME_REVEAL/P5_V2_REVEAL_CORE_RESULTS.tsv')
fac=rows('13_MPT_OUTCOME_REVEAL/P5_V2_REVEAL_FACTOR_RESULTS.tsv')
rep=rows('13_MPT_OUTCOME_REVEAL/P5_V2_REVEAL_REPLICATE_RESULTS.tsv')
gene=rows('13_MPT_OUTCOME_REVEAL/P5_V2_GENE_LEVEL_COMPANION_RESULTS.tsv')
osta=rows('14_OSTA_SPECIFICITY/P5_V2_OSTA_SPECIFICITY_CORE_RESULTS.tsv')
spatial=rows(EVO/'03_FULL_EVOLUTION_DATA/P5_STAGE1B_FIRST_BRIDGE_MODEL_RESULTS.tsv')
sp_summary=source(EVO/'03_FULL_EVOLUTION_DATA/P5_STAGE1B_FIRST_FROZEN_BRIDGE_REVEAL_SUMMARY.txt').read_text()
audit=source('13_MPT_OUTCOME_REVEAL/P5_V2_IMPLEMENTATION_AUDIT_v1.log').read_text()
def extract(key):return re.search(r'^'+re.escape(key)+r' = (.+)$',audit,re.M).group(1)
shared=[]
for prefix,label in [('SHARED_PRIMARY','Downstream +30..+90'),('SHARED_NEGATIVE','Earlier -90..-30')]:
    ci=extract(prefix+'_CI').strip('[]').split(', ')
    shared.append(dict(analysis=label,N=extract('SHARED_PRIMARY_NEGATIVE_TMD'),genes=extract('SHARED_PRIMARY_NEGATIVE_GENES'),beta=extract(prefix+'_BETA'),ci_low=ci[0],ci_high=ci[1],p=extract(prefix+'_P'),se=''))
repro=rows('15_PUBLICATION_VALIDATION/01_AUTHOR_REPRODUCTION/P5_V2_AUTHOR_REPRO_PROFILE_RESULTS.tsv')[:9]
matched=rows('15_PUBLICATION_VALIDATION/01_AUTHOR_REPRODUCTION/P5_V2_AUTHOR_REPRO_MATCHED_VALUES.tsv.gz')
attr=rows('15_PUBLICATION_VALIDATION/02_ATTRITION_AUDIT/P5_V2_INCLUDED_EXCLUDED_DESCRIPTIVES.tsv')
flow=rows('15_PUBLICATION_VALIDATION/02_ATTRITION_AUDIT/P5_V2_ATTRITION_FLOW.tsv')
multi=rows('15_PUBLICATION_VALIDATION/02_ATTRITION_AUDIT/INTERNAL_TMD_MULTIPLICITY_DISTRIBUTIONS.tsv')
bconst=rows('15_PUBLICATION_VALIDATION/03_B_NORMALIZATION/P5_V2_B_CONSTANTS.tsv')
binv=rows('15_PUBLICATION_VALIDATION/03_B_NORMALIZATION/P5_V2_B_INVARIANCE_RESULTS.tsv')

# Tables copy archived character strings; no numerical recomputation.
st1=[]
for data,path,group in [(core,'13_MPT_OUTCOME_REVEAL/P5_V2_REVEAL_CORE_RESULTS.tsv','MPT core'),(fac,'13_MPT_OUTCOME_REVEAL/P5_V2_REVEAL_FACTOR_RESULTS.tsv','MPT factor'),(rep,'13_MPT_OUTCOME_REVEAL/P5_V2_REVEAL_REPLICATE_RESULTS.tsv','MPT profile'),(gene,'13_MPT_OUTCOME_REVEAL/P5_V2_GENE_LEVEL_COMPANION_RESULTS.tsv','Gene companion'),(osta,'14_OSTA_SPECIFICITY/P5_V2_OSTA_SPECIFICITY_CORE_RESULTS.tsv','Machinery comparison'),(shared,'13_MPT_OUTCOME_REVEAL/P5_V2_IMPLEMENTATION_AUDIT_v1.log','Shared temporal audit')]:
    for idx,r in enumerate(data):
        label=r.get('analysis',r.get('factor',r.get('profile','')))
        cls='PRIMARY' if group=='MPT core' and label=='PRIMARY' else 'PRESPECIFIED_SECONDARY'
        timing='Specified before MPT primary reveal'
        note='TMD observations; gene-clustered inference; score standardized over 5538 bridge TMDs.'
        if group=='MPT profile':cls='SUPPORTING';note+=' Descriptive constituent profile; not an independent replication.'
        if group=='Gene companion':cls='SUPPORTING';timing='Broad role pre-reveal; detailed implementation post-primary reveal';note='Gene observations; HC3 inference; gene means standardized over 1400 bridge genes. Blank Holm P for aggregate is not zero.'
        if group=='Machinery comparison':cls='POST_REVEAL_SECONDARY';timing='Broad OST-A role pre-MPT reveal; detailed shared-population contrast specified after MPT adjudication';note+=' Individual P values are not a machinery-difference test.'
        if group=='Shared temporal audit':cls='SUPPORTING';timing='Earlier corridor prespecified; shared primary estimate added in post-reveal descriptive audit';note+=' No between-window test. SE absent from this audit is left blank; printed precision preserved.'
        st1.append(dict(table_group=group,evidence_class=cls,timing=timing,source_file=path,source_location=('section 7' if group=='Shared temporal audit' else 'data row '+str(idx+1)),**r,interpretation_note=note))
tsv('supplement/ST1.tsv',st1)
st2=[]
for idx,r in enumerate(spatial):
    within=r['model']=='WITHIN_PROTEIN_D_MEAN'
    note='Historical orthogonal evidence. Location coefficients secondary to the joint test; raw P must be read with Holm P.'
    if within:note+=' Specification–implementation discrepancy: executed model demeaned outcome and score within gene–location groups; specification contained different covariate-handling language, including covariate demeaning. No harmonized reanalysis was performed. SUPPORTING only.'
    st2.append(dict(evidence_class='HISTORICAL_SUPPORTING' if within or r['model']!='PRIMARY_D_MEAN' else 'HISTORICAL_ORTHOGONAL',**r,source_file='04_STAGE1_EVOLUTION/03_FULL_EVOLUTION_DATA/P5_STAGE1B_FIRST_BRIDGE_MODEL_RESULTS.tsv',source_location='data row '+str(idx+1),interpretation_note=note))
tsv('supplement/ST2.tsv',st2)
joint_specs=[('Primary joint','5.3786305','3','0.14608079','2158','487'),('Compartment difference','1.7434982','2','0.41821939','2158','487'),('Adjusted joint','5.6800092','3','0.12825984','2158','487'),('Within-protein joint','5.9226784','3','0.11543329','1625','265'),('Replicate 1 joint','','','0.32692243','',''),('Replicate 2 joint','','','0.089193294','',''),('Pseudo component joint','5.4427739','3','0.14210172','',''),('Real component joint','2.0319536','3','0.56580117','','')]
for row in joint_specs:
    for val in row[1:]:
        if val:assert val in sp_summary
tsv('supplement/ST2_joint_tests.tsv',[dict(analysis=a,Wald=w,df=d,p=p,N_stacked=n,gene_clusters=g,evidence_class='HISTORICAL_SUPPORTING' if a not in ['Primary joint','Compartment difference'] else 'HISTORICAL_ORTHOGONAL',note='Within-protein result has the specification–implementation discrepancy disclosed in ST2. Blank cells are not supplied in this summary and are not inferred.',source_file='04_STAGE1_EVOLUTION/03_FULL_EVOLUTION_DATA/P5_STAGE1B_FIRST_FROZEN_BRIDGE_REVEAL_SUMMARY.txt') for a,w,d,p,n,g in joint_specs])
for name,data in [('S1_profile_statistics',repro),('S2_attrition_flow',flow),('S2_attributes',attr),('S2_multiplicity',multi),('S3_B_constants',bconst),('S3_invariance',binv),('Figure3_temporal',shared)]:tsv('source_data/'+name+'.tsv',data)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
blue='#245878';orange='#A75923';gray='#52606D'
inventory=[]
def finish(fig,name,title):
    fig.suptitle(title,x=.04,ha='left',fontsize=15,fontweight='bold')
    for ext in ['pdf','svg','png']:
        fig.savefig(OUT/'figures'/f'{name}.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig);inventory.append(name)
def panel(ax,label,title):ax.set_title(label+'  '+title,loc='left',fontweight='bold',fontsize=11,pad=13)
def forest(ax,data,labels,letter,title,keys=('beta','ci_low','ci_high'),pkey='p',xlim=(-.13,.09)):
    for i,(r,label) in enumerate(zip(data,labels)):
        beta,lo,hi=[float(r[k]) for k in keys]
        ax.plot([lo,hi],[i,i],color=blue,lw=1.8);ax.plot(beta,i,'o',color=blue,ms=5)
        ax.text(1.02,i,'P='+format(float(r[pkey]),'.6g') if r.get(pkey) else '',transform=ax.get_yaxis_transform(),fontsize=9,va='center')
    ax.axvline(0,color='#9CA3AF',ls='--',lw=.8);ax.set_yticks(range(len(data)),labels);ax.set_ylim(len(data)-.4,-.6);ax.set_xlim(*xlim);ax.grid(axis='x',alpha=.15)
    ax.set_xlabel('Score coefficient (95% CI)');panel(ax,letter,title)
def textpanel(ax,letter,title,lines):
    ax.axis('off');panel(ax,letter,title)
    for i,(head,body) in enumerate(lines):
        y=.91-i*(.85/max(len(lines),1));ax.text(0,y,head,va='top',fontsize=12,color=blue,fontweight='bold');ax.text(0,y-.08,body,va='top',fontsize=10,linespacing=1.5)

fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
textpanel(axs[0,0],'A','Evolutionary predictor · definition',[('Internal TMD C-terminal end → +30..+90','61 codons; retained evolutionary annotation boundaries'),('Conserved amino-acid positions: A, G, P, T, V','Synonymous third-base matches versus same-gene,\nsame-amino-acid expectation'),('Z = (observed − expected) / √variance','Median across ≥6 valid species; positive = excess conservation.\nNot a codon-optimality or translation-rate measurement.')])
textpanel(axs[0,1],'B','Evolutionary coverage · descriptive',[('8,420 internal TMDs','6,523 with complete corridor geometry'),('6,511 scored TMDs','1,656 genes; frozen before functional outcome inspection')])
textpanel(axs[1,0],'C','Protein-locked correspondence · supporting',[('Exact full-protein identity','Unambiguous reciprocal, order-preserving TMD correspondence'),('5,538 scored TMDs / 1,400 genes','Corresponding boundaries need not be identical.\nEach annotation retains its own +30..+90 anchor.')])
textpanel(axs[1,1],'D','Measurable population · descriptive',[('5,538 mapped TMDs / 1,400 genes','All mapped functional windows pass geometry'),('2,796 measurable TMDs / 736 genes','Positive paired-input baseline in all six profiles;\n50.4875% TMD retention'),('2,742 excluded TMDs / 664 genes','Measurement availability; no enrichment threshold selection')])
finish(fig,'Figure_1','Evolutionary predictor and measurable population')

fig=plt.figure(figsize=(12,10));gs=fig.add_gridspec(3,1,height_ratios=[3,3,5],left=.25,right=.83,top=.91,bottom=.07,hspace=.65)
forest(fig.add_subplot(gs[0]),core[:3],['Primary · PRIMARY','Adjusted · SECONDARY','Within-protein · SECONDARY'],'A','MPT hierarchy; gene-clustered inference')
forest(fig.add_subplot(gs[1]),fac,[r['factor'] for r in fac],'B','Factor estimates · SECONDARY; Holm P',pkey='holm_p')
forest(fig.add_subplot(gs[2]),rep,[r['profile'].replace('_',' ') for r in rep],'C','Constituent profiles · SUPPORTING, not independent replication')
finish(fig,'Figure_2','MPT-associated ribosome enrichment: estimates and uncertainty')

fig=plt.figure(figsize=(12,10));gs=fig.add_gridspec(3,1,height_ratios=[4,2,3],left=.26,right=.83,top=.91,bottom=.07,hspace=.75)
forest(fig.add_subplot(gs[0]),gene,['Combined MPT','TMCO1','CCDC47','Nicalin'],'A','Gene-level companion · SUPPORTING; 652 genes')
forest(fig.add_subplot(gs[1]),shared,[r['analysis'] for r in shared],'B','Temporal audit · SUPPORTING; 2,455 TMDs / 642 genes')
forest(fig.add_subplot(gs[2]),osta,['MPT shared','OST-A shared','MPT minus OST-A'],'C','Machinery comparison · POST-REVEAL SECONDARY; 2,575 / 683')
finish(fig,'Figure_3','Gene-level, temporal and machinery boundaries')

fig=plt.figure(figsize=(12,8));gs=fig.add_gridspec(2,1,height_ratios=[1,1.25],left=.2,right=.83,top=.89,bottom=.07,hspace=.55)
ax=fig.add_subplot(gs[0]);forest(ax,spatial[:3],['ER','Mitochondria','Cytosol'],'A','Historical spatial analysis; 2,158 stacked observations / 487 genes',keys=('beta_per_1SD_evolution_score','CI_2.5','CI_97.5'),pkey='Holm_p',xlim=(-.1,.2));ax.set_xlabel('Spatial score coefficient (95% CI); Holm P at right')
ax.text(0,1.02,'Joint P=0.14608079; compartment-difference P=0.41821939',transform=ax.transAxes,fontsize=10)
textpanel(fig.add_subplot(gs[1]),'B','Integrated hierarchy · interpretive synthesis',[('Observed evidence','Negative constituent MPT directions; prespecified secondary\ncovariate-adjusted association; supporting measurement fidelity'),('Not established','Primary MPT and within-protein association; gene-level and\nbroader spatial generalization; temporal or machinery specificity'),('Interpretation retained','A weaker or broader evolutionary-functional link remains possible.\nNo causal or elongation-rate inference; no equivalence or biological absence.')])
finish(fig,'Figure_4','Historical orthogonal evidence and limits of interpretation')

fig,axs=plt.subplots(3,3,figsize=(12,11),layout='constrained')
for ax,r in zip(axs.flat,repro):
    pts=[m for m in matched if m['profile']==r['profile'] and m['group']==r['group']]
    ax.scatter([float(m['author_enrichment']) for m in pts],[float(m['reconstructed_enrichment']) for m in pts],s=2,alpha=.18,color=blue,rasterized=True)
    ax.axline((0,0),slope=1,color=gray,lw=.8,ls='--')
    ax.set_title(r['profile']+'  N='+r['matched_N'],loc='left',fontsize=11)
    ax.text(.04,.96,'Pearson r='+format(float(r['Pearson_r']),'.6f')+'\nSaved slope='+format(float(r['OLS_slope_reconstruction_on_author']),'.6f')+'\nIntercept='+format(float(r['OLS_intercept']),'.6f'),transform=ax.transAxes,va='top',fontsize=8)
    ax.set_xlabel('Author log2 enrichment');ax.set_ylabel('Reconstructed log2 enrichment')
finish(fig,'Figure_S1','S1  Author reconstruction fidelity · SUPPORTING; dashed line = identity')
# Saved calibration summaries are supplied without any fit or new quantile.
fig,axs=plt.subplots(1,2,figsize=(12,6),gridspec_kw={'width_ratios':[1.25,1]});fig.subplots_adjust(left=.22,right=.97,top=.83,bottom=.15,wspace=.6)
labels=['Evolutionary score','Amino-acid identity','GC3','CpG density','TMD length','Hydropathy','TMD index','Downstream length','Downstream lumenal']
axs[0].barh(range(len(attr)),[float(r['SMD_included_minus_excluded']) for r in attr],color=blue);axs[0].set_yticks(range(len(attr)),labels);axs[0].invert_yaxis();axs[0].axvline(0,color=gray,lw=.8);axs[0].set_xlabel('Saved included − excluded SMD');panel(axs[0],'A','Attribute differences · DESCRIPTIVE')
for key,label,c in [('all_prelocked_genes_prelocked_TMDs_per_gene','All mapped genes',gray),('formal_valid_TMDs_per_contributing_gene','Contributing genes',blue)]:
    rr=[r for r in multi if r['distribution']==key];axs[1].plot([int(r['internal_TMDs_per_gene']) for r in rr],[float(r['fraction']) for r in rr],marker='o',ms=4,label=label,color=c)
axs[1].set_xlabel('Internal TMDs per gene');axs[1].set_ylabel('Saved fraction of genes');axs[1].legend(frameon=False,fontsize=8);panel(axs[1],'B','Multiplicity · DESCRIPTIVE')
finish(fig,'Figure_S2','S2  Measurement availability: 5,538 → 2,796 TMDs; no attrition tests')

fig,axs=plt.subplots(1,2,figsize=(12,6));fig.subplots_adjust(left=.15,right=.96,top=.84,bottom=.2,wspace=.6)
axs[0].barh(range(9),[float(r['expected_log2_shift']) for r in bconst],color=blue);axs[0].set_yticks(range(9),[r['profile'].replace('_',' ') for r in bconst]);axs[0].invert_yaxis();axs[0].set_xlabel('Saved additive log2 shift');panel(axs[0],'A','Profile-wide B constants · SUPPORTING')
rr=[r for r in binv if r['kind']=='association']
for field,label,m in [('beta_delta','Score coefficient','o'),('se_delta','Clustered SE','s'),('p_delta','P value','^')]:
    axs[1].plot([float(r[field]) for r in rr],range(len(rr)),m,label=label,ms=5)
axs[1].set_yticks(range(len(rr)),[r['analysis'].replace('_',' ') for r in rr],fontsize=8);axs[1].invert_yaxis();axs[1].axvline(0,color=gray,lw=.8);axs[1].set_xlabel('Saved full-B minus B-free difference');axs[1].legend(frameon=False,fontsize=8,loc='lower left');panel(axs[1],'B','Archived invariance audit · SUPPORTING')
fig.text(.15,.07,'Primary intercept shifts by the common constant; primary score inference is preserved.\nSaved combined MPT shift: 0.50673694561000371; nonconstant residual: 8.8817841970012523e-16.',fontsize=9)
finish(fig,'Figure_S3','S3  B-normalization technical validation; no biological reanalysis')

tsv('INPUT_PROVENANCE_SHA256.tsv',[dict(source_file=p,sha256=h) for p,h in sources.items()])
print('Rendered archived values only; no manuscript edits or scientific fits.')

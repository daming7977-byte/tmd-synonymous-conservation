from pathlib import Path
import sys,subprocess,os,json
out=Path(__file__).resolve().parent
sys.path.insert(0,str(out/'dependencies'))
os.environ['PYTHONPATH']=str(out/'dependencies');os.environ['PYTHONDONTWRITEBYTECODE']='1'
import pysam
q=out/'counting_fixture';q.mkdir(exist_ok=True)
(q/'fixture.gtf').write_text('chrT\ttest\tCDS\t101\t120\t.\t+\t0\tgene_id "A"; transcript_id "A1";\nchrT\ttest\tCDS\t111\t130\t.\t+\t0\tgene_id "B"; transcript_id "B1";\nchrT\ttest\tCDS\t101\t120\t.\t-\t0\tgene_id "C"; transcript_id "C1";\nchrT\ttest\tCDS\t101\t120\t.\t+\t0\tgene_id "A"; transcript_id "A2";\n')
# A-only, A+B ambiguous, reverse-strand C, low MAPQ excluded,
# and a splice spanning the CDS entirely in the skipped N block.
records=[('a',95,0,50,[(0,10)]),('ab',110,0,50,[(0,10)]),('c',105,16,50,[(0,10)]),('low',100,0,0,[(0,10)]),('splice',90,0,50,[(0,5),(3,40),(0,5)])]
with pysam.AlignmentFile(str(q/'fixture.bam'),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chrT','LN':1000}]}) as f:
 for name,start,flag,mq,cigar in sorted(records,key=lambda x:x[1]):
  r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='A'*10;r.flag=flag;r.reference_id=0;r.reference_start=start;r.mapping_quality=mq;r.cigartuples=cigar;r.set_tag('NH',1);f.write(r)
cmd=[sys.executable,'-m','HTSeq.scripts.count','--format=bam','--order=pos','--stranded=yes','--type=CDS','--idattr=gene_id','--mode=union','--nonunique=all','--minaqual=10',str(q/'fixture.bam'),str(q/'fixture.gtf')]
p=subprocess.run(cmd,capture_output=True,text=True,check=True)
actual={a:int(b) for a,b in [l.split('\t') for l in p.stdout.splitlines()]}
assert all(actual[g]==v for g,v in {'A':2,'B':1,'C':1,'__too_low_aQual':1,'__no_feature':1,'__ambiguous':1}.items()),actual
(q/'verification.json').write_text(json.dumps({'command':cmd,'counts':actual,'verified':'same strand; nonunique all; no duplicate isoform multiplication; MAPQ default; skipped CIGAR region excluded'},indent=2)+'\n')
print('HTSeq 2.0.3 counting fixture PASS',actual)

import gzip,os,sys
DB="/sessions/ecstatic-optimistic-planck/mnt/Susanth/gingicharm/src/pgmlst/data"
ROOT="/sessions/ecstatic-optimistic-planck/mnt/Susanth"
LOCI=["ftsQ","gpdxJ","hagB","mcmA","pepO","pga","recA"]
G={"1439":"Pg_clinical_isolates/1439.fna.gz",
"KCOM_2799":"Pg_Genomes/clinical_isolates/KCOM_2799__GCF_002753955.1.fasta",
"Kyudai-3":"Pg_clinical_isolates/Kyudai-3.fna","Kyudai-4":"Pg_clinical_isolates/Kyudai-4.fna",
"SJD12":"Pg_clinical_isolates/SJD12.fna.gz","WW2866":"Pg_clinical_isolates/WW2866.fna",
"JKG10":"Pg_clinical_isolates/JKG10.fna","LyEC01":"Pg_Genomes/clinical_isolates/LyEC01__GCF_030144345.1.fasta",
"NBRC_115148":"Pg_clinical_isolates/NBRC_115148.fna","WW2096":"Pg_clinical_isolates/WW2096.fna"}
ST={"1439":208,"KCOM_2799":209,"Kyudai-3":210,"Kyudai-4":211,"SJD12":212,
    "WW2866":213,"JKG10":217,"LyEC01":218,"NBRC_115148":219,"WW2096":220}
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def readfas(p):
    op=gzip.open if p.endswith(".gz") else open
    d={};name=None;buf=[]
    for line in op(p,"rt"):
        line=line.strip()
        if line.startswith(">"):
            if name: d[name]="".join(buf)
            name=line[1:];buf=[]
        elif line: buf.append(line)
    if name: d[name]="".join(buf)
    return d
al={}
for loc in LOCI:
    al[loc]={k.split("_")[1]:v.upper() for k,v in readfas(os.path.join(DB,loc+".fas")).items()}
prof={}
for l in open(os.path.join(DB,"profiles.tsv")).readlines()[1:]:
    f=l.rstrip("\n").split("\t"); prof[int(f[0])]=tuple(f[1:8])
bad=0
for strain,st in sorted(ST.items(), key=lambda kv: kv[1]):
    seq="N".join(readfas(os.path.join(ROOT,G[strain])).values()).upper()
    seqrc=rc(seq.replace("N","A")); seq2=seq
    found=[]
    for loc in LOCI:
        hits=sorted((n for n,s in al[loc].items() if s in seq2 or rc(s) in seq2), key=int)
        found.append(hits)
    expect=prof[st]
    okloci=[e in h for e,h in zip(expect,found)]
    multi=[f"{l}={'/'.join(h)}" for l,h,e in zip(LOCI,found,expect) if len(h)>1]
    status="OK " if all(okloci) else "BAD"
    if not all(okloci): bad+=1
    print(f"{status} ST{st:<4}{strain:<14}expected {'-'.join(expect)}  recovered {'-'.join('/'.join(h) or '-' for h in found)}"
          + (f"   [multi-copy: {', '.join(multi)}]" if multi else ""))
print(f"\n{len(ST)-bad}/{len(ST)} profiles re-derived from the assembly by exact allele match")
sys.exit(1 if bad else 0)

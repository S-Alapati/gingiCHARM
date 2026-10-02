import os,gzip,itertools,collections,csv
DB="/sessions/ecstatic-optimistic-planck/mnt/Susanth/gingicharm/src/pgmlst/data"
ROOT="/sessions/ecstatic-optimistic-planck/mnt/Susanth"
LOCI=["ftsQ","gpdxJ","hagB","mcmA","pepO","pga","recA"]
EXPLEN={"ftsQ":420,"gpdxJ":380,"hagB":380,"mcmA":420,"pepO":400,"pga":400,"recA":400}
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def readfas(p):
    op=gzip.open if p.endswith(".gz") else open
    d={};n=None;b=[]
    for l in op(p,"rt"):
        l=l.strip()
        if l.startswith(">"):
            if n: d[n]="".join(b)
            n=l[1:];b=[]
        elif l: b.append(l)
    if n: d[n]="".join(b)
    return d
al={loc:{k.split("_")[1]:v.upper() for k,v in readfas(os.path.join(DB,loc+".fas")).items()} for loc in LOCI}
SEED=24
seeds=collections.defaultdict(list)
for loc in LOCI:
    for n,s in al[loc].items(): seeds[s[:SEED]].append((loc,n,s))
prof={}
for l in open(os.path.join(DB,"profiles.tsv")).readlines()[1:]:
    f=l.rstrip("\n").split("\t"); prof[tuple(f[1:8])]=f[0]
OLD=os.path.join(ROOT,"pgtoolkit/docs/SupplementaryTableS3_MLST_calls.tsv")
old=[l.rstrip("\n").split("\t") for l in open(OLD) if l.strip()]
oldrows=old[1:]
idx={}
for d in ["Pg_clinical_isolates","Pg_Genomes/clinical_isolates","Pg_Genomes/reference_strains","Pg_Genomes","Pg_Unified_DB","Pg_CoreGenome_DB"]:
    dd=os.path.join(ROOT,d)
    if not os.path.isdir(dd): continue
    for f in os.listdir(dd):
        if f.endswith((".fna",".fna.gz",".fasta",".fasta.gz",".fa")): idx.setdefault(f,os.path.join(dd,f))
OURS=set(str(x) for x in list(range(178,200))+[208,209,210,211,212,213,217,218,219,220])
out=[]; st_status=collections.Counter(); notfound=[]
for r in oldrows:
    strain,acc=r[0],r[1]
    path=None
    for fn,p in idx.items():
        base=fn.split("__")[0]
        for suf in (".fna.gz",".fasta.gz",".fasta",".fna",".fa"): base=base.replace(suf,"")
        if base==strain or (acc and acc in fn): path=p; break
    if not path:
        notfound.append((strain,acc)); out.append(r); continue
    fwd="N".join(readfas(path).values()).upper(); both=fwd+"NN"+rc(fwd)
    hit=collections.defaultdict(set)
    for seed,entries in seeds.items():
        start=0
        while True:
            i=both.find(seed,start)
            if i<0: break
            for loc,n,s in entries:
                if both[i:i+len(s)]==s: hit[loc].add(n)
            start=i+1
    cands=[sorted(hit.get(loc,()),key=int) for loc in LOCI]
    multi=[l for l,h in zip(LOCI,cands) if len(h)>1]
    missing=[l for l,h in zip(LOCI,cands) if not h]
    pstr="-".join("/".join(h) if h else "-" for h in cands)
    # distinguish undeposited allele from locus not recovered, using the previous call
    undep=[x for x in (r[5].split(",") if len(r)>5 and r[5] else []) if x]
    notrec=[x for x in (r[6].split(",") if len(r)>6 and r[6] else []) if x]
    undep=[l for l in undep if l in missing]; notrec=[l for l in notrec if l in missing]
    unexplained=[l for l in missing if l not in undep and l not in notrec]
    undep+=unexplained
    if missing:
        st=""; status=("locus not recovered end to end" if notrec and not undep
                       else "undeposited allele")
    else:
        hits={prof[c] for c in itertools.product(*cands) if c in prof}
        if len(hits)==1:
            n=sorted(hits,key=int)[0]; st="ST"+n
            status="defined from this work" if n in OURS else "deposited"
        elif len(hits)>1:
            st=""; status="unresolved multi-copy locus"
        else:
            st=""; status="all alleles deposited, type not yet defined"
    st_status[status]+=1
    out.append([strain,acc,st,status,pstr,",".join(sorted(set(undep))),
                ",".join(sorted(set(notrec))),",".join(multi)])
with open(OLD,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(old[0]); w.writerows(out)
print("genomes not located:",notfound)
print("\nstatus breakdown over",len(out),"assemblies:")
for k,v in st_status.most_common(): print(f"  {v:>4}  {k}")
defined=[r[2][2:] for r in out if len(r)>2 and r[2].startswith("ST")]
print(f"\n{len(defined)} of {len(out)} assemblies carry a defined ST, spanning {len(set(defined))} distinct types")
print(f"{len(set(defined)&OURS)} of those types were defined by this work")
print("multi-copy hagB seen in", sum(1 for r in out if len(r)>7 and "hagB" in r[7]), "assemblies")
print("\nall alleles deposited, type not yet defined:")
for r in out:
    if len(r)>3 and r[3].startswith("all alleles"): print("   ",r[0],r[4])

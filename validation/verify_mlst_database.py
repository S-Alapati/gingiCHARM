import os,re,sys,collections
DB="/sessions/ecstatic-optimistic-planck/mnt/Susanth/gingicharm/src/pgmlst/data"
SRC="/sessions/ecstatic-optimistic-planck/mnt/Susanth/gingicharm/data/pubmlst_source"
LOCI=["ftsQ","gpdxJ","hagB","mcmA","pepO","pga","recA"]
EXPLEN={"ftsQ":420,"gpdxJ":380,"hagB":380,"mcmA":420,"pepO":400,"pga":400,"recA":400}
fails=[]
def ck(cond,msg):
    print(("  PASS  " if cond else "  FAIL  ")+msg)
    if not cond: fails.append(msg)

def readfas(p):
    d={};name=None;buf=[]
    for line in open(p):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if name: d[name]="".join(buf)
            name=line[1:].strip();buf=[]
        elif line.strip(): buf.append(line.strip())
    if name: d[name]="".join(buf)
    return d

print("\n1. ALLELE FILES")
alleles={}; seqtoid={}
for loc in LOCI:
    d=readfas(os.path.join(DB,loc+".fas"))
    nums=sorted(int(k.split("_")[1]) for k in d)
    ck(nums==list(range(1,len(nums)+1)), f"{loc}: {len(nums)} alleles numbered 1..{len(nums)} with no gaps")
    bad=[k for k,v in d.items() if len(v)!=EXPLEN[loc]]
    ck(not bad, f"{loc}: every allele is {EXPLEN[loc]} bp"+(f" (bad: {bad})" if bad else ""))
    nonacgt=[k for k,v in d.items() if not re.fullmatch(r"[ACGTacgt]+",v)]
    ck(not nonacgt, f"{loc}: ACGT only"+(f" ({nonacgt})" if nonacgt else ""))
    rev=collections.defaultdict(list)
    for k,v in d.items(): rev[v.upper()].append(k)
    dup={s:v for s,v in rev.items() if len(v)>1}
    ck(not dup, f"{loc}: no two alleles share a sequence"+(f" ({list(dup.values())})" if dup else ""))
    alleles[loc]={k.split("_")[1]:v.upper() for k,v in d.items()}
    seqtoid[loc]={v.upper():k.split("_")[1] for k,v in d.items()}
tot=sum(len(v) for v in alleles.values())
ck(tot==241, f"total alleles = {tot} (expect 241)")

print("\n2. PROFILE TABLE")
rows=[l.rstrip("\n").split("\t") for l in open(os.path.join(DB,"profiles.tsv")) if l.strip()]
hdr=rows[0]; prof={}
ck(hdr[:8]==["ST"]+LOCI, f"header = {hdr[:8]}")
for r in rows[1:]: prof[int(r[0])]=tuple(r[1:8])
ck(len(prof)==210, f"{len(prof)} profiles (expect 210)")
ck(len(set(prof.values()))==len(prof), "every profile is unique")
missingal=[(st,l,a) for st,p in prof.items() for l,a in zip(LOCI,p) if a not in alleles[l]]
ck(not missingal, "every allele referenced by a profile exists"+(f" ({missingal[:5]})" if missingal else ""))
gaps=[x for x in range(1,221) if x not in prof]
ck(gaps==[201,202,203,204,205,206,207,214,215,216],
   f"only the ten unconfirmed PRJDB2925 STs are absent: {gaps}")
new={208:"1439",209:"KCOM_2799",210:"Kyudai-3",211:"Kyudai-4",212:"SJD12",
     213:"WW2866",217:"JKG10",218:"LyEC01",219:"NBRC_115148",220:"WW2096"}
ck(all(st in prof for st in new), "all ten new STs present")

print("\n3. NEW PROFILES AGAINST THE SUBMITTED FILE")
sub={l.split("\t")[0]:tuple(l.rstrip("\n").split("\t")[1:8])
     for l in open("/sessions/ecstatic-optimistic-planck/mnt/Susanth/PubMLST_SUBMISSION_FINAL/update4_combined_profiles/profiles_upload.tsv").readlines()[1:] if l.strip()}
for st,strain in sorted(new.items()):
    ck(prof[st]==sub[strain], f"ST{st} = {'-'.join(prof[st])} matches submitted {strain}")

print("\n4. PUBMLST 221-ISOLATE VALIDATION")
iso={}
for l in open(os.path.join(SRC,"pubMLST_isolateIDs.txt")).readlines()[1:]:
    f=l.rstrip("\n").split("\t")
    if len(f)>11: iso[f[0]]=dict(name=f[1], profile=tuple(f[4:11]), st=f[11].strip().rstrip("="))
ck(len(iso)==221, f"{len(iso)} isolate records loaded")
seqs=collections.defaultdict(dict)
for h,s in readfas(os.path.join(SRC,"pg_mlst_sequences.txt")).items():
    m=re.match(r"(\d+):\S+\s+\S+\s+(\w+)",h)
    if m: seqs[m.group(1)][m.group(2)]=s.upper().rstrip("=")  # WW3039 carries a trailing = in the PubMLST export
ck(len(seqs)==221, f"{len(seqs)} isolates have sequence data")
st_by_prof={v:k for k,v in prof.items()}
ok=0; mism=[]
for iid,rec in iso.items():
    called=[]
    for loc in LOCI:
        s=seqs.get(iid,{}).get(loc)
        called.append(seqtoid[loc].get(s) if s else None)
    if None in called: mism.append((iid,rec["name"],"locus not matched",called)); continue
    if tuple(called)!=rec["profile"]: mism.append((iid,rec["name"],"profile differs",tuple(called),rec["profile"])); continue
    got=st_by_prof.get(tuple(called))
    if rec["st"] and str(got)!=rec["st"]: mism.append((iid,rec["name"],"ST differs",got,rec["st"])); continue
    ok+=1
ck(ok==221, f"{ok}/221 isolates re-typed to their deposited allele profile and ST")
for m in mism[:10]: print("      ",m)

print("\n5. CROSS-PACKAGE SYNC")
PC="/sessions/ecstatic-optimistic-planck/mnt/Susanth/periocharm/src/periocore/data/mlst/P_gingivalis"
for loc in LOCI+["profiles"]:
    f=(loc+".fas") if loc!="profiles" else "profiles.tsv"
    a=open(os.path.join(DB,f)).read(); b=open(os.path.join(PC,f)).read()
    ck(a==b, f"{f} identical in gingiCHARM and perioCHARM")

print("\n"+("ALL CHECKS PASSED" if not fails else f"{len(fails)} CHECK(S) FAILED"))
sys.exit(1 if fails else 0)

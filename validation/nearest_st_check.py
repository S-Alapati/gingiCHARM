import csv,os
LOCI=["ftsQ","gpdxJ","hagB","mcmA","pepO","pga","recA"]
DB="/sessions/ecstatic-optimistic-planck/mnt/Susanth/gingicharm/src/pgmlst/data"
def profiles():
    with open(os.path.join(DB,"profiles.tsv")) as fh:
        return list(csv.DictReader(fh,delimiter="\t"))
PROF=profiles()
def nearest(cands,max_ties=4):
    scored=[]
    for pr in PROF:
        hit=[l for l in LOCI if cands.get(l) and pr[l] in cands[l]]
        scored.append((len(hit),int(pr["ST"]),pr,[l for l in LOCI if l not in hit]))
    scored.sort(key=lambda x:(-x[0],x[1]))
    top=scored[0][0]; tied=[s for s in scored if s[0]==top]; f=tied[0]
    return {"ST":f[2]["ST"],"loci_matched":top,"mismatches":len(LOCI)-top,
            "mismatch_loci":f[3],"ties":[t[2]["ST"] for t in tied[1:max_ties+1]],
            "n_tied":len(tied)}
# ---- validation against the PubMLST query the user ran ----
ww={"ftsQ":[],"gpdxJ":["35"],"hagB":["1"],"mcmA":["3"],"pepO":["8"],"pga":["7"],"recA":["8"]}
n=nearest(ww)
print("WW2885, against PubMLST's own answer")
print(f"  PubMLST : ST-111 and 2 others, 3 mismatches, 4/7 (57.1 %)")
print(f"  pgmlst  : ST{n['ST']} and {n['n_tied']-1} others, {n['mismatches']} mismatches, "
      f"{n['loci_matched']}/7 ({100*n['loci_matched']/7:.1f} %)")
ok = n["ST"]=="111" and n["n_tied"]==3 and n["mismatches"]==3 and n["loci_matched"]==4
print("  ->", "MATCHES" if ok else "DIFFERS")
print("  mismatching loci:", ", ".join(n["mismatch_loci"]))
print("  ties:", ", ".join("ST"+t for t in n["ties"]))
# ---- exact-match sanity: a defined ST must score 7/7 on itself ----
bad=0
for pr in PROF[:40]:
    c={l:[pr[l]] for l in LOCI}
    nn=nearest(c)
    if nn["loci_matched"]!=7 or nn["ST"]!=pr["ST"]: bad+=1
print(f"\nself-consistency on 40 defined profiles: {40-bad}/40 return themselves at 7/7")

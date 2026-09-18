# PubMLST source files

Retrieved from the *P. gingivalis* PubMLST databases on 17 September 2026 and
kept for provenance. The package itself reads `src/pgmlst/data/`, not these.

| File | Contents |
|---|---|
| `PG_MLST_Profiles.txt` | 200 seven-locus profiles, ST1 to ST200 |
| `pubMLST_isolateIDs.txt` | 221 isolate records with their allele designations |
| `pg_mlst_sequences.txt` | concatenated seven-locus sequence for each of the 221 isolates |

`pg_mlst_sequences.txt` is what makes independent validation possible: running
the caller over it reproduces the curated designations and the sequence type for
219 of the 221 records, with no disagreements. The two exceptions are isolate
196 (COM30Pg1.1), typed with alleles that are assigned but not yet public, and
isolate 221 (WW3039), whose exported sequences each carry a trailing `=`
character that has to be stripped before they will resolve.

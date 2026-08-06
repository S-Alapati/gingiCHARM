# pgtoolkit_web

A small Flask web GUI for [pgtoolkit](../pgtoolkit/). Upload one or more
*P. gingivalis* genome assemblies, tick the analyses you want
(virulence-factor screening, fimA / mfa-operon / rag / K-antigen typing,
functional categories), and view the results in a browser — or download
them as TSV.

The backend is the same `pgvfdb`, `fimtyping`, `mfatyping`, `ragtyping`,
`kantigen` and `pgfunc` modules the pgtoolkit package ships. The web app
calls them, collects the results, and renders them as tables.

## Running it

1. Install pgtoolkit first — from the sibling `pgtoolkit/` folder:

   ```
   cd ../pgtoolkit
   bash install.sh
   ```

   This sets up a virtualenv, installs pgtoolkit, pulls in NCBI BLAST+, and
   bundles the reference database. The web app reuses the same venv.

2. Launch the web app:

   ```
   bash start.sh
   ```

   It will install Flask if it isn't there yet, then start the server on
   <http://127.0.0.1:5057>. Open that URL in any browser.

3. From the landing page click **Start an analysis →**, drop one or more
   FASTAs into the upload box, tick the modules you want, and hit **Run
   analysis**. For a single genome you'll get an inline per-module results
   page; for multiple, a batch summary table that puts every genome on one
   row.

## What you can do from the GUI

- **Upload multiple files at once** — drag a folder of clinical isolates in
  and analyse them in one go.
- **Pick analyses individually** — six checkboxes, one per module. Default
  is everything selected.
- **Tune the BLAST thresholds** — minimum identity, minimum coverage and
  thread count are exposed; defaults are 60% identity / 50% coverage / 4
  threads (same as pgtoolkit's CLI).
- **Download as TSV** — one button for all results, one button per module.
  The TSV is long-format (`genome / module / headline / field / value`) so
  it merges cleanly across runs.

## Layout

```
pgtoolkit_web/
├── app.py                    Flask backend
├── start.sh                  one-shot launcher
├── requirements.txt          Flask + Werkzeug
├── templates/
│   ├── welcome.html          landing page
│   ├── index.html            upload + module selection
│   ├── results.html          single-genome results
│   ├── batch_results.html    multi-genome results
│   └── about.html
├── uploads/                  uploaded FASTAs (auto-created)
└── results/                  saved run payloads as JSON (auto-created)
```

## Re-opening past runs

Every run gets a job ID (e.g. `20260602_133012_a1b2c3`) and a JSON file in
`results/`. Past runs are reachable at
`http://127.0.0.1:5057/result/<job_id>` — no re-uploading or
re-analysis needed.

## Notes and limitations

- The web app runs every selected module sequentially on each uploaded
  genome. A full six-module run on a complete assembly takes about 30–60
  seconds on a modern laptop.
- Uploads are capped at 200 MB total per submission. Increase
  `MAX_CONTENT_LENGTH` in `app.py` if you need bigger.
- The server binds to `127.0.0.1` only. To expose it on the network,
  change the `host` argument at the bottom of `app.py`.
- `PORT=8080 bash start.sh` runs it on a different port.

## When something goes wrong

- *"Cannot find pgtoolkit"* on startup — pgtoolkit isn't installed.
  Either run `bash install.sh` inside `pgtoolkit/` first, or set
  `PGTOOLKIT_SRC` to the path of `pgtoolkit/src` and re-run.
- *"BLAST executable not found"* — pgtoolkit installs `ncbi-blast+`
  automatically, but if the install was skipped or PATH isn't set,
  install BLAST+ separately or re-run the pgtoolkit installer.
- *"upload over the 200 MB limit"* — split the batch, or raise the limit
  in `app.py`.

"""
gingiCHARM — Flask GUI front-end for the pgtoolkit Python package.

Upload one or more P. gingivalis genome assemblies, pick which typing
modules to run (virulence factors, fimA, mfa operon, rag, K-antigen,
functional categories), and view the results inline or download as TSV.

This is a thin wrapper around the `pgtoolkit` Python package. The package
must be installed (run `bash install.sh` from the pgtoolkit/ folder) or
reachable via the PGTOOLKIT_SRC env var.
"""
import csv
import io
import json
import os
import sys
import time
import uuid
from datetime import datetime

# Locate pgtoolkit. If it's installed, this just works; otherwise look one
# directory up from this file for the repository's src/ tree.
def _bootstrap_pgtoolkit():
    try:
        import pgcore  # noqa: F401
        return
    except ImportError:
        pass
    candidates = []
    here = os.path.dirname(os.path.abspath(__file__))
    if os.environ.get("PGTOOLKIT_SRC"):
        candidates.append(os.environ["PGTOOLKIT_SRC"])
    candidates += [
        os.path.join(here, "..", "src"),          # repo layout: web/ beside src/
        os.path.join(here, "..", "pgtoolkit", "src"),
        os.path.join(here, "pgtoolkit", "src"),
    ]
    for c in candidates:
        if os.path.isdir(os.path.join(c, "pgcore")):
            sys.path.insert(0, os.path.abspath(c))
            return
    raise RuntimeError(
        "Cannot find the gingiCHARM package. Run 'bash install.sh' from the "
        "repository root, or set PGTOOLKIT_SRC to the path of src/.")


_bootstrap_pgtoolkit()

try:
    from flask import (Flask, render_template, request, send_file, abort,
                       redirect, url_for)
    from werkzeug.utils import secure_filename
except ImportError:
    sys.stderr.write(
        "\n"
        "  pgtoolkit_web needs Flask, which isn't installed in this Python.\n"
        "\n"
        "  Install it with one of:\n"
        "     pip install flask                 # if pip is on PATH\n"
        "     python3 -m pip install flask      # always works\n"
        "     conda install -c conda-forge flask\n"
        "\n"
        "  Then re-run:  python3 app.py     (or:  bash start.sh)\n"
        "\n")
    sys.exit(1)

import pgvfdb
import fimtyping
import mfatyping
import ragtyping
import kantigen
import pgfunc

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024  # 200 MB / upload
HERE = os.path.dirname(os.path.abspath(__file__))
app.config["UPLOAD_FOLDER"] = os.path.join(HERE, "uploads")
app.config["RESULTS_FOLDER"] = os.path.join(HERE, "results")
os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)
os.makedirs(app.config["RESULTS_FOLDER"], exist_ok=True)

# Map module key -> (label, callable, description)
MODULES = {
    "pgvfdb":    ("Virulence factors",       pgvfdb.analyze,
                  "Screen the genome against 154 curated P. gingivalis virulence-factor families."),
    "fimtyping": ("fimA typing",             fimtyping.analyze,
                  "Assign the FimA fimbrial type (I–V) from the full-length fimA gene."),
    "mfatyping": ("mfa operon typing",       mfatyping.analyze,
                  "Type the whole minor-fimbriae operon (mfa1/mfa2/mfa3/mfa4/mfa5) following Watanabe 2021 / Nagano 2023."),
    "ragtyping": ("rag typing",              ragtyping.analyze,
                  "Type the rag pathogenicity-island locus (rag-1, rag-3, rag-4) via ragA + ragB."),
    "kantigen":  ("K-antigen typing",        kantigen.analyze,
                  "Assign the capsular K-antigen serotype from the cps locus (K1, K4, K6, K-)."),
    "pgfunc":    ("Functional categories",   pgfunc.analyze,
                  "Place each gene into one of ten functional categories."),
}


def _run_jobs(genome_path, selected_keys, min_identity, min_coverage, threads):
    """Run the selected modules against one genome. Returns a list of dicts."""
    jobs = []
    for key in selected_keys:
        if key not in MODULES: continue
        label, fn, _ = MODULES[key]
        t0 = time.time()
        try:
            res = fn(genome_path, min_identity=min_identity,
                     min_coverage=min_coverage, threads=threads)
            jobs.append({
                "key": key, "label": label, "ok": True,
                "headline": getattr(res, "headline", ""),
                "summary":  getattr(res, "summary", {}),
                "columns":  list(res.columns),
                "records":  list(res.records),
                "elapsed":  round(time.time() - t0, 2),
            })
        except Exception as e:
            jobs.append({
                "key": key, "label": label, "ok": False,
                "error": f"{type(e).__name__}: {e}",
                "elapsed": round(time.time() - t0, 2),
            })
    return jobs


def _save_run(job_id, payload):
    path = os.path.join(app.config["RESULTS_FOLDER"], f"{job_id}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    return path


def _load_run(job_id):
    path = os.path.join(app.config["RESULTS_FOLDER"], f"{job_id}.json")
    if not os.path.exists(path): return None
    with open(path) as fh:
        return json.load(fh)


# Routes -------------------------------------------------------------

@app.route("/")
def welcome():
    return render_template("welcome.html")


@app.route("/analyze", methods=["GET", "POST"])
def analyze():
    if request.method == "GET":
        return render_template("index.html", modules=MODULES)

    files = [f for f in request.files.getlist("genome") if f and f.filename]
    if not files:
        return render_template("index.html", modules=MODULES,
                               error="Pick at least one genome FASTA to upload.")
    selected = request.form.getlist("modules")
    if not selected:
        return render_template("index.html", modules=MODULES,
                               error="Select at least one analysis to run.")

    try:
        min_identity = float(request.form.get("min_identity", "60"))
        min_coverage = float(request.form.get("min_coverage", "50"))
        threads      = int(request.form.get("threads", "4"))
    except ValueError:
        return render_template("index.html", modules=MODULES,
                               error="min identity / coverage / threads must be numbers.")

    job_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:6]
    genomes = []
    for f in files:
        safe = secure_filename(f.filename)
        dest = os.path.join(app.config["UPLOAD_FOLDER"], f"{job_id}__{safe}")
        f.save(dest)
        genomes.append({"name": safe, "path": dest})

    runs = []
    for g in genomes:
        jobs = _run_jobs(g["path"], selected, min_identity, min_coverage, threads)
        runs.append({"genome": g["name"], "jobs": jobs})

    payload = {
        "job_id": job_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "min_identity": min_identity,
        "min_coverage": min_coverage,
        "threads": threads,
        "selected": selected,
        "runs": runs,
    }
    _save_run(job_id, payload)

    template = "batch_results.html" if len(runs) > 1 else "results.html"
    return render_template(template, payload=payload, modules=MODULES)


@app.route("/result/<job_id>")
def result(job_id):
    payload = _load_run(job_id)
    if payload is None: abort(404)
    template = "batch_results.html" if len(payload["runs"]) > 1 else "results.html"
    return render_template(template, payload=payload, modules=MODULES)


@app.route("/download/<job_id>/<scope>.tsv")
def download_tsv(job_id, scope):
    """Download analysis results as TSV.
    scope = "all" (everything), or "<module_key>" (one module).
    """
    payload = _load_run(job_id)
    if payload is None: abort(404)
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t")
    w.writerow(["genome", "module", "headline", "field", "value"])
    for run in payload["runs"]:
        for j in run["jobs"]:
            if scope != "all" and j["key"] != scope: continue
            if not j["ok"]:
                w.writerow([run["genome"], j["label"], "(error)", "error", j["error"]])
                continue
            head = j.get("headline", "")
            for k, v in (j.get("summary") or {}).items():
                w.writerow([run["genome"], j["label"], head, k, v])
            for rec in j["records"]:
                for col in j["columns"]:
                    w.writerow([run["genome"], j["label"], head, col,
                                rec.get(col, "")])
    data = buf.getvalue().encode("utf-8")
    return send_file(io.BytesIO(data), mimetype="text/tab-separated-values",
                     as_attachment=True,
                     download_name=f"pgtoolkit_{job_id}_{scope}.tsv")


@app.route("/about")
def about():
    return render_template("about.html")


@app.errorhandler(413)
def too_large(e):
    return render_template("index.html", modules=MODULES,
                           error="That upload is over the 200 MB limit. "
                                 "Try fewer / smaller files."), 413


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5057"))
    print(f"\n  gingiCHARM running at  http://127.0.0.1:{port}\n")
    app.run(host="127.0.0.1", port=port, debug=False)

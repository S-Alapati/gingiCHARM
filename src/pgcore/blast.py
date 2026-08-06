"""Locating and running NCBI BLAST+."""
import glob
import os
import shutil
import subprocess

_HOME = os.path.expanduser("~/.pgtoolkit/blast")


def _candidate_dirs():
    dirs = []
    env = os.environ.get("PGTOOLKIT_BLAST_BIN")
    if env:
        dirs.append(env)
    dirs += sorted(glob.glob(os.path.join(_HOME, "*", "bin")))
    dirs += [os.path.join(_HOME, "bin"), _HOME]
    return dirs


def find_blast(program="blastn"):
    """Return the path to a BLAST+ executable, or raise a helpful error."""
    for d in _candidate_dirs():
        p = os.path.join(d, program)
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    p = shutil.which(program)
    if p:
        return p
    raise RuntimeError(
        "Could not find '%s' (NCBI BLAST+).\n"
        "Run the pgtoolkit installer (install.sh) to download BLAST+, or set\n"
        "the PGTOOLKIT_BLAST_BIN environment variable to its bin directory."
        % program)


def blast_available():
    """True if BLAST+ can be found."""
    try:
        find_blast("blastn")
        return True
    except RuntimeError:
        return False


_FIELDS = ("qseqid sseqid pident length mismatch gapopen qstart qend "
           "sstart send evalue bitscore slen qlen")
_INT = ("length", "slen", "qlen", "qstart", "qend", "sstart", "send",
        "mismatch", "gapopen")


def run_blast(program, query, db, *, evalue=1e-5, threads=4,
              max_target_seqs=50, extra=None):
    """Run a BLAST search and return a list of HSP dicts.

    Each dict carries the standard outfmt-6 fields plus 'scov' (subject
    coverage = aligned length / subject length, as a percentage).
    """
    exe = find_blast(program)
    cmd = [exe, "-query", query, "-db", db, "-outfmt", "6 " + _FIELDS,
           "-evalue", str(evalue), "-num_threads", str(threads),
           "-max_target_seqs", str(max_target_seqs)]
    if extra:
        cmd += list(extra)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("%s failed:\n%s" % (program, proc.stderr.strip()))
    cols = _FIELDS.split()
    hits = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != len(cols):
            continue
        h = dict(zip(cols, parts))
        for k in ("pident", "evalue", "bitscore"):
            h[k] = float(h[k])
        for k in _INT:
            h[k] = int(h[k])
        h["scov"] = round(100.0 * h["length"] / h["slen"], 1) if h["slen"] else 0.0
        hits.append(h)
    return hits

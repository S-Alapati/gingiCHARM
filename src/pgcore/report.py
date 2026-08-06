"""Result container and shared command-line runner."""
import argparse
import csv
import sys


class Result:
    """A uniform result object returned by every pgtoolkit analysis."""

    def __init__(self, title, columns, records, summary=None, headline=None):
        self.title = title
        self.columns = list(columns)
        self.records = list(records)        # list of dict
        self.summary = dict(summary or {})  # dict of headline metrics
        self.headline = headline            # one-line key finding

    def __len__(self):
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    def __repr__(self):
        lines = ["<%s: %d rows>" % (self.title, len(self.records))]
        if self.headline:
            lines.append("  " + self.headline)
        for k, v in self.summary.items():
            lines.append("  %s: %s" % (k, v))
        return "\n".join(lines)

    def to_tsv(self, path):
        """Write all records to a tab-separated file."""
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=self.columns, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            for r in self.records:
                w.writerow(r)
        return path

    def to_dataframe(self):
        """Return the records as a pandas DataFrame (pandas must be installed)."""
        import pandas as pd
        return pd.DataFrame(self.records, columns=self.columns)

    def print_table(self, limit=None, stream=sys.stdout):
        rows = self.records if limit is None else self.records[:limit]
        if not rows:
            print("  (no rows)", file=stream)
            return
        w = {c: len(c) for c in self.columns}
        for r in rows:
            for c in self.columns:
                w[c] = max(w[c], len(str(r.get(c, ""))))
        print("  ".join(c.ljust(w[c]) for c in self.columns), file=stream)
        print("  ".join("-" * w[c] for c in self.columns), file=stream)
        for r in rows:
            print("  ".join(str(r.get(c, "")).ljust(w[c]) for c in self.columns),
                  file=stream)
        if limit is not None and len(self.records) > limit:
            print("  ... %d more rows (use -o to write the full table)"
                  % (len(self.records) - limit), file=stream)


def run_cli(prog, description, analyze_func, *, min_identity=80.0,
            min_coverage=70.0, argv=None):
    """Shared argparse-based command-line entry point for every module."""
    p = argparse.ArgumentParser(prog=prog, description=description)
    p.add_argument("query",
                   help="FASTA file: genome assembly, contigs, CDS or proteins")
    p.add_argument("-o", "--output", metavar="TSV",
                   help="write the full result table to this file")
    p.add_argument("--min-identity", type=float, default=min_identity,
                   help="minimum %% identity (default: %s)" % min_identity)
    p.add_argument("--min-coverage", type=float, default=min_coverage,
                   help="minimum %% coverage (default: %s)" % min_coverage)
    p.add_argument("--threads", type=int, default=4, help="BLAST threads")
    a = p.parse_args(argv)
    try:
        res = analyze_func(a.query, min_identity=a.min_identity,
                           min_coverage=a.min_coverage, threads=a.threads)
    except Exception as e:                                  # noqa: BLE001
        print("error: %s" % e, file=sys.stderr)
        return 1
    print("# " + res.title)
    if res.headline:
        print("# " + res.headline)
    print()
    res.print_table(limit=50)
    if res.summary:
        print()
        for k, v in res.summary.items():
            print("  %s: %s" % (k, v))
    if a.output:
        res.to_tsv(a.output)
        print("\nfull table -> %s" % a.output)
    return 0

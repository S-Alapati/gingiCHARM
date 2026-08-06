#!/usr/bin/env bash
# pgtoolkit installer - sets up the package, its modules and all dependencies
# (NCBI BLAST+). Works on macOS and Linux.
#
#   bash install.sh
#
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"

PY=""
for c in python3 python; do
    if command -v "$c" >/dev/null 2>&1; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
    echo "error: Python 3 is required but was not found on PATH." >&2
    echo "Install Python 3.8+ from https://www.python.org/downloads/ and retry." >&2
    exit 1
fi

echo "pgtoolkit installer - using $("$PY" --version 2>&1)"
exec "$PY" "$HERE/install.py" "$@"

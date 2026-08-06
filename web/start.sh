#!/bin/bash
# Launch the gingiCHARM web GUI.
# Run this from the repository's web/ folder, after installing the package
# with `bash install.sh` at the repository root.

set -e
HERE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT="$( cd "$HERE/.." && pwd )"
cd "$HERE"

PYTHON="${PYTHON:-python3}"

# Reuse the package venv created by install.sh at the repo root.
if [ -d "$ROOT/.venv" ]; then
  # shellcheck disable=SC1091
  source "$ROOT/.venv/bin/activate"
  PYTHON=python
  echo "Using virtualenv at $ROOT/.venv"
fi

# Install Flask if it isn't already
if ! $PYTHON -c "import flask" 2>/dev/null; then
  echo "Installing Flask into $($PYTHON -c 'import sys; print(sys.prefix)') ..."
  if ! $PYTHON -m pip install --quiet -r requirements.txt 2>/dev/null; then
    $PYTHON -m pip install --quiet --break-system-packages -r requirements.txt || {
      echo ""
      echo "Could not auto-install Flask. Install it manually with:"
      echo "    $PYTHON -m pip install flask"
      echo "or, if you use conda:"
      echo "    conda install -c conda-forge flask"
      exit 1
    }
  fi
fi

# Make the package importable; fall back to the repo's src/ tree.
if ! $PYTHON -c "import pgcore" 2>/dev/null; then
  if [ -d "$ROOT/src/pgcore" ]; then
    export PGTOOLKIT_SRC="$ROOT/src"
    echo "Using package source at $PGTOOLKIT_SRC"
  else
    echo "ERROR: the gingiCHARM package is not installed and src/ was not"
    echo "       found at $ROOT/src. Run 'bash install.sh' from the"
    echo "       repository root first."
    exit 1
  fi
fi

# Locate NCBI BLAST+.
if ! command -v blastn >/dev/null 2>&1 && [ -z "$PGTOOLKIT_BLAST_BIN" ]; then
  for cand in \
      "$ROOT/.venv/bin" \
      "$ROOT/blast/bin" \
      "$ROOT/ncbi-blast/bin" \
      /usr/local/ncbi/blast/bin \
      /opt/homebrew/opt/blast/bin \
      /opt/homebrew/bin; do
    if [ -x "$cand/blastn" ]; then
      export PGTOOLKIT_BLAST_BIN="$cand"
      echo "Found BLAST+ at $PGTOOLKIT_BLAST_BIN"
      break
    fi
  done
fi

PORT="${PORT:-5057}"
echo ""
echo "  gingiCHARM"
echo "  Open  http://127.0.0.1:$PORT  in your browser"
echo "  Press Ctrl+C to stop."
echo ""

PORT=$PORT $PYTHON app.py

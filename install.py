#!/usr/bin/env python3
"""
pgtoolkit installer.

Sets up everything needed to run the toolkit:
  1. creates an isolated Python virtual environment (~/.pgtoolkit/venv)
  2. installs the pgtoolkit package (the six analysis modules) into it
  3. downloads NCBI BLAST+ (the only external dependency) into ~/.pgtoolkit/blast
  4. links the command-line tools into ~/.local/bin

Works on macOS and Linux. Run with:  python3 install.py
"""
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
HOME = os.path.expanduser("~/.pgtoolkit")
VENV = os.path.join(HOME, "venv")
BLAST_DIR = os.path.join(HOME, "blast")
BIN_LINK = os.path.expanduser("~/.local/bin")
BLAST_VERSION = "2.17.0"
COMMANDS = ["pgvfdb", "fimtyping", "mfatyping", "ragtyping", "kantigen",
            "pgfunc", "pgtoolkit-info"]


def log(msg):
    print("[pgtoolkit-install] " + msg, flush=True)


def blast_platform_tags():
    """Candidate NCBI BLAST+ platform tags for this machine, best first."""
    sysname, mach = platform.system(), platform.machine().lower()
    if sysname == "Darwin":
        if mach in ("arm64", "aarch64"):
            return ["aarch64-macosx", "x64-macosx"]
        return ["x64-macosx"]
    if sysname == "Linux":
        if mach in ("aarch64", "arm64"):
            return ["aarch64-linux", "x64-linux"]
        return ["x64-linux"]
    raise RuntimeError("Unsupported OS: %s (install BLAST+ manually)" % sysname)


def have_blast():
    """True if a BLAST+ blastn is already available to the toolkit."""
    import glob
    dirs = glob.glob(os.path.join(BLAST_DIR, "*", "bin")) + [BLAST_DIR]
    env = os.environ.get("PGTOOLKIT_BLAST_BIN")
    if env:
        dirs.insert(0, env)
    for d in dirs:
        if os.path.isfile(os.path.join(d, "blastn")):
            return True
    return shutil.which("blastn") is not None


def download_blast():
    if have_blast():
        log("BLAST+ already available - skipping download.")
        return
    os.makedirs(BLAST_DIR, exist_ok=True)
    base = "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/%s/" % BLAST_VERSION
    last_err = None
    for tag in blast_platform_tags():
        fname = "ncbi-blast-%s+-%s.tar.gz" % (BLAST_VERSION, tag)
        url = base + fname
        log("downloading NCBI BLAST+ (%s) ... this is ~250 MB" % tag)
        try:
            tmp = os.path.join(tempfile.gettempdir(), fname)
            urllib.request.urlretrieve(url, tmp)
            log("extracting BLAST+ ...")
            with tarfile.open(tmp) as tf:
                tf.extractall(BLAST_DIR)
            os.remove(tmp)
            log("BLAST+ installed in %s" % BLAST_DIR)
            return
        except Exception as e:                                  # noqa: BLE001
            last_err = e
            log("  could not fetch %s (%s)" % (tag, e))
    raise RuntimeError(
        "Failed to download BLAST+ automatically (%s).\n"
        "Install it manually from "
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ and set the\n"
        "PGTOOLKIT_BLAST_BIN environment variable to its bin/ directory."
        % last_err)


def make_venv_and_install():
    if not os.path.isdir(VENV):
        log("creating virtual environment: %s" % VENV)
        subprocess.run([sys.executable, "-m", "venv", VENV], check=True)
    pip = os.path.join(VENV, "bin", "pip")
    log("upgrading pip ...")
    subprocess.run([pip, "install", "-q", "--upgrade", "pip"], check=True)
    log("installing pgtoolkit (six analysis modules) ...")
    subprocess.run([pip, "install", "-q", PKG_DIR], check=True)
    log("pgtoolkit package installed.")


def link_commands():
    os.makedirs(BIN_LINK, exist_ok=True)
    linked = []
    for cmd in COMMANDS:
        src = os.path.join(VENV, "bin", cmd)
        dst = os.path.join(BIN_LINK, cmd)
        if os.path.isfile(src):
            if os.path.islink(dst) or os.path.exists(dst):
                try:
                    os.remove(dst)
                except OSError:
                    pass
            try:
                os.symlink(src, dst)
                linked.append(cmd)
            except OSError:
                pass
    if linked:
        log("linked commands into %s : %s" % (BIN_LINK, ", ".join(linked)))


def main():
    log("installing into %s" % HOME)
    os.makedirs(HOME, exist_ok=True)
    make_venv_and_install()
    try:
        download_blast()
    except RuntimeError as e:
        log("WARNING: " + str(e))
    link_commands()

    # verify
    info = os.path.join(VENV, "bin", "pgtoolkit-info")
    print()
    if os.path.isfile(info):
        env = dict(os.environ)
        import glob
        for d in glob.glob(os.path.join(BLAST_DIR, "*", "bin")):
            env["PGTOOLKIT_BLAST_BIN"] = d
        subprocess.run([info], env=env)

    print()
    log("DONE.")
    print("""
Use the toolkit in either of two ways:

  1. Command line - the tools are in %s
        pgvfdb     my_assembly.fna
        fimtyping  my_assembly.fna
        mfatyping  my_assembly.fna
        ragtyping  my_assembly.fna
        kantigen   my_assembly.fna
        pgfunc     my_assembly.fna
     (add that folder to your PATH if the commands are not found:
        export PATH="$HOME/.local/bin:$PATH")

  2. Python API - activate the environment first
        source %s/bin/activate
        python -c "import pgvfdb; print(pgvfdb.analyze('my_assembly.fna'))"
""" % (BIN_LINK, VENV))


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        log("ERROR: a command failed (%s)" % e)
        sys.exit(1)
    except Exception as e:                                      # noqa: BLE001
        log("ERROR: %s" % e)
        sys.exit(1)

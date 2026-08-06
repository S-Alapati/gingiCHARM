"""FASTA parsing, translation and nucleotide/protein detection."""

_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODON = {a + b + c: _AA[i] for i, (a, b, c) in
          enumerate((x, y, z) for x in _BASES for y in _BASES for z in _BASES)}


def read_fasta(path):
    """Return {id: sequence}; id is the first whitespace-delimited token."""
    seqs, h, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(buf)
                tok = line[1:].split()
                h, buf = (tok[0] if tok else line[1:]), []
            elif h is not None:
                buf.append(line.strip())
        if h is not None:
            seqs[h] = "".join(buf)
    return seqs


def translate(nt):
    """Translate a nucleotide sequence in frame 1 (standard genetic code)."""
    nt = nt.upper().replace("U", "T")
    return "".join(_CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def guess_type(seqs):
    """Guess 'nt' or 'aa' from a sequence string or a {id: seq} dict."""
    if isinstance(seqs, dict):
        sample = "".join(list(seqs.values())[:50])
    else:
        sample = seqs
    sample = sample.upper()[:20000]
    if not sample:
        return "nt"
    acgtn = sum(c in "ACGTNU" for c in sample)
    return "nt" if acgtn / len(sample) > 0.9 else "aa"

### Boas Pucker ###
### bpucker@cebitec.uni-bielefelde ###
#!/usr/bin/env python3
# v0.3

__usage__ = """
python transeq.py --in <INPUT_FILE|INPUT_DIR> --out <OUTPUT_FILE|OUTPUT_DIR> [--internal-stop-to-x]

Notes:
- Unknown/ambiguous codons (e.g., containing N) translate to 'X'.
- If --internal-stop-to-x is provided, internal '*' are converted to 'X' (a final '*' is kept).
"""

import sys, os, glob

# --- helpers --- #

def load_multiple_fasta_file(fasta_file):
    """Load all sequences from a (possibly wrapped) FASTA file into a dict."""
    content = {}
    header = None
    seq_chunks = []
    with open(fasta_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(">"):
                if header is not None:
                    content[header] = "".join(seq_chunks)
                header = line[1:].strip().split()[0]  # trim after first whitespace
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            content[header] = "".join(seq_chunks)
    return content

def translate(seq, genetic_code, unknown="X", internal_stop_to_x=False, keep_terminal_stop=True):
    """Translate DNA to protein.
       - Unknown/ambiguous codons -> 'X'
       - Stop codons per table -> '*'
       - Optionally convert internal '*' to 'X' (keeps final '*' if present).
    """
    seq = seq.upper().replace("U", "T")
    pep = []
    n = len(seq) // 3
    for i in range(n):
        codon = seq[i*3:i*3+3]
        aa = genetic_code.get(codon)
        if aa is None:
            aa = unknown
        pep.append(aa)
    pep = "".join(pep)

    if internal_stop_to_x and "*" in pep:
        if keep_terminal_stop and pep.endswith("*"):
            pep = pep[:-1].replace("*", "X") + "*"
        else:
            pep = pep.replace("*", "X")
    return pep

def translate_file(in_fa, out_fa, genetic_code, internal_stop_to_x=False):
    seqs = load_multiple_fasta_file(in_fa)
    internal_stop_count = 0
    with open(out_fa, "w") as out:
        for header, nt in seqs.items():
            pep = translate(nt, genetic_code, unknown="X",
                            internal_stop_to_x=internal_stop_to_x,
                            keep_terminal_stop=True)
            if "*" in pep[:-1]:
                internal_stop_count += 1
            out.write(f">{header}\n{pep}\n")
    print(f"[INFO] {os.path.basename(in_fa)} -> {os.path.basename(out_fa)} | "
          f"seqs: {len(seqs)} | internal-stops (pre-fix): {internal_stop_count}")

def gather_inputs(spec):
    """Return list of input files. If directory, grab common FASTA extensions."""
    if os.path.isdir(spec):
        exts = ("*.fa", "*.fasta", "*.fna", "*.fas")
        files = []
        for ext in exts:
            files.extend(glob.glob(os.path.join(spec, ext)))
        files.sort()
        return files
    else:
        return [spec]

def make_output_path(out_spec, in_file):
    """Decide output path: if out_spec is a dir or endswith '/', write there; else treat as file path."""
    if out_spec.endswith(os.sep) or os.path.isdir(out_spec):
        os.makedirs(out_spec, exist_ok=True)
        base = os.path.basename(in_file)
        root = os.path.splitext(base)[0]
        return os.path.join(out_spec, f"{root}.pep.fa")
    else:
        # Single input -> exact output file path
        # Multi-input with a single-file out_spec is not supported
        return out_spec

def main(argv):
    if '--in' not in argv or '--out' not in argv:
        sys.exit(__usage__)

    input_spec = argv[argv.index('--in')+1]
    output_spec = argv[argv.index('--out')+1]
    internal_stop_to_x = ('--internal-stop-to-x' in argv)

    # Standard code table (nuclear, 1)
    genetic_code = {
        'CTT':'L','ATG':'M','AAG':'K','AAA':'K','ATC':'I','AAC':'N','ATA':'I','AGG':'R',
        'CCT':'P','ACT':'T','AGC':'S','ACA':'T','AGA':'R','CAT':'H','AAT':'N','ATT':'I',
        'CTG':'L','CTA':'L','CTC':'L','CAC':'H','ACG':'T','CCG':'P','AGT':'S','CAG':'Q',
        'CAA':'Q','CCC':'P','TAG':'*','TAT':'Y','GGT':'G','TGT':'C','CGA':'R','CCA':'P',
        'TCT':'S','GAT':'D','CGG':'R','TTT':'F','TGC':'C','GGG':'G','TGA':'*','GGA':'G',
        'TGG':'W','GGC':'G','TAC':'Y','GAG':'E','TCG':'S','TTA':'L','GAC':'D','TCC':'S',
        'GAA':'E','TCA':'S','GCA':'A','GTA':'V','GCC':'A','GTC':'V','GCG':'A','GTG':'V',
        'TTC':'F','GTT':'V','GCT':'A','ACC':'T','TTG':'L','CGT':'R','TAA':'*','CGC':'R'
    }

    inputs = gather_inputs(input_spec)

    # If multiple inputs but output_spec is a file path, abort to avoid overwriting
    if len(inputs) > 1 and (not output_spec.endswith(os.sep) and not os.path.isdir(output_spec)):
        sys.exit("[ERROR] For multiple input files, --out must be a directory or end with '/'.")
    if output_spec.endswith(os.sep):
        os.makedirs(output_spec, exist_ok=True)

    for in_file in inputs:
        out_file = make_output_path(output_spec, in_file)
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        translate_file(in_file, out_file, genetic_code, internal_stop_to_x=internal_stop_to_x)

if __name__ == '__main__':
    main(sys.argv)

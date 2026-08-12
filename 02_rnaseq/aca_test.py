import gzip, csv, statistics

seqs, name, buf = {}, None, []
with gzip.open('ref/ecoli_cds.fna.gz','rt') as f:
    for line in f:
        if line.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)

def aca_per_kb(s):
    s = s.upper()
    n = sum(1 for i in range(len(s)-2) if s[i:i+3] == 'ACA')
    return n / len(s) * 1000

rows = []
with open('results/deseq2_results.csv') as f:
    for x in csv.DictReader(f):
        tid, padj = x[''], x['padj']
        if tid in seqs and padj not in ('NA',''):
            s = seqs[tid]
            if len(s) < 150: continue
            rows.append((float(x['log2FoldChange']), aca_per_kb(s), float(padj)))

down = [d for lfc,d,p in rows if p < 0.05 and lfc < -1]
up   = [d for lfc,d,p in rows if p < 0.05 and lfc >  1]
ns   = [d for lfc,d,p in rows if p >= 0.05]

for label, grp in [('DOWN (sig, LFC<-1)', down), ('UP (sig, LFC>1)', up), ('Not significant', ns)]:
    if grp:
        print(f"{label:22s} n={len(grp):5d}  median ACA/kb = {statistics.median(grp):.2f}  mean = {statistics.mean(grp):.2f}")

# --- confounder check: does ACA density track expression level? ---
import math
rows2 = []
with open('results/deseq2_results.csv') as f:
    for x in csv.DictReader(f):
        tid = x['']
        if tid in seqs and x['padj'] not in ('NA',''):
            s = seqs[tid]
            if len(s) < 150: continue
            rows2.append((float(x['baseMean']), aca_per_kb(s), float(x['log2FoldChange'])))

rows2.sort(key=lambda r: r[0])
n = len(rows2); q = n // 4
print("\n--- ACA density by expression quartile ---")
for i, label in enumerate(['Q1 lowest','Q2','Q3','Q4 highest']):
    chunk = rows2[i*q:(i+1)*q] if i < 3 else rows2[3*q:]
    print(f"{label:12s} n={len(chunk):5d}  median ACA/kb = {statistics.median([r[1] for r in chunk]):.2f}  median LFC = {statistics.median([r[2] for r in chunk]):+.2f}")

# Mann-Whitney U, DOWN vs NS
try:
    from scipy.stats import mannwhitneyu
    u, p = mannwhitneyu(down, ns, alternative='greater')
    print(f"\nMann-Whitney (DOWN > NS): p = {p:.3g}")
except ImportError:
    print("\n(scipy not installed - run: mamba install -y scipy)")

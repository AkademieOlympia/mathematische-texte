from collections import Counter
from math import log2
import re

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def parse_genbank(path):
    text = open(path, 'r', encoding='utf-8').read()
    locus = re.search(r'^LOCUS\s+(\S+)', text, re.M)
    accession = re.search(r'^ACCESSION\s+(\S+)', text, re.M)
    version = re.search(r'^VERSION\s+(\S+)', text, re.M)
    definition = re.search(r'^DEFINITION\s+(.+?)(?=^ACCESSION)', text, re.M | re.S)
    organism = re.search(r'^\s+ORGANISM\s+(.+)$', text, re.M)

    origin = re.search(r'ORIGIN(.*)//', text, re.S)
    if origin is None:
        raise ValueError('Keine ORIGIN-Sequenz gefunden.')
    seq = ''.join(re.findall(r'[ACGTUacgtu]+', origin.group(1))).upper().replace('U', 'T')

    features = []
    lines = text.splitlines()
    in_features = False
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('FEATURES'):
            in_features = True
            i += 1
            continue
        if in_features and line.startswith('ORIGIN'):
            break
        if in_features and re.match(r'^\s{5}CDS\s+\d+\.\.\d+', line):
            m = re.match(r'^\s{5}CDS\s+(\d+)\.\.(\d+)', line)
            start, end = int(m.group(1)), int(m.group(2))
            quals = {}
            i += 1
            current_key = None
            while i < len(lines):
                qline = lines[i]
                if re.match(r'^\s{5}\S', qline) or qline.startswith('ORIGIN'):
                    i -= 1
                    break
                qm = re.match(r'^\s{21}/([^=]+)="?(.*?)"?\s*$', qline)
                if qm:
                    current_key = qm.group(1)
                    value = qm.group(2)
                    quals[current_key] = value
                else:
                    cont = qline.strip()
                    if current_key and cont:
                        quals[current_key] = quals.get(current_key, '') + cont.strip('"')
                i += 1
            features.append({'start': start, 'end': end, 'qualifiers': quals})
        i += 1

    return {
        'locus': locus.group(1) if locus else None,
        'accession': accession.group(1) if accession else None,
        'version': version.group(1) if version else None,
        'definition': ' '.join(definition.group(1).split()) if definition else None,
        'organism': organism.group(1) if organism else None,
        'sequence': seq,
        'features': features,
    }


def gc_content(seq):
    return 100.0 * sum(1 for b in seq if b in 'GC') / len(seq) if seq else 0.0


def shannon_entropy(seq):
    counts = Counter(seq)
    n = len(seq)
    return -sum((c/n) * log2(c/n) for c in counts.values()) if n else 0.0


def kmer_counts(seq, k):
    return Counter(seq[i:i+k] for i in range(len(seq)-k+1))


def codons(seq):
    usable = len(seq) - (len(seq) % 3)
    return [seq[i:i+3] for i in range(0, usable, 3)]


def translate(seq):
    return ''.join(CODON_TABLE.get(c, 'X') for c in codons(seq))


def base_counts_by_codon_position(seq):
    pos = [Counter(), Counter(), Counter()]
    usable = len(seq) - (len(seq) % 3)
    for i in range(usable):
        pos[i % 3][seq[i]] += 1
    return pos


def transition_matrix(seq):
    bases = 'ACGT'
    M = matrix(ZZ, 4, 4, 0)
    idx = {b: i for i, b in enumerate(bases)}
    for a, b in zip(seq, seq[1:]):
        if a in idx and b in idx:
            M[idx[a], idx[b]] += 1
    return bases, M


def print_top(counter, n=10, title='Top-Einträge'):
    print(title)
    for key, value in counter.most_common(n):
        print('  {:>8} : {}'.format(key, value))


def analyze_feature(seq, feature):
    start, end = feature['start'], feature['end']
    sub = seq[start-1:end]
    aa = translate(sub)
    print('\nCDS {}..{}   Gen: {}'.format(start, end, feature['qualifiers'].get('gene', '?')))
    print('  Länge nt: {}'.format(len(sub)))
    print('  Länge Codons: {}'.format(len(codons(sub))))
    print('  GC-Gehalt: {:.2f}%'.format(gc_content(sub)))
    print('  Protein: {}'.format(aa))
    annotated = feature['qualifiers'].get('translation', '').replace(' ', '')
    if annotated:
        ok = aa.rstrip('*') == annotated
        print('  Annotation stimmt mit Translation überein: {}'.format(ok))


def main(path='DNA.gb'):
    data = parse_genbank(path)
    seq = data['sequence']

    print('=== GenBank / SageMath-Analyse ===')
    print('Locus      :', data['locus'])
    print('Accession  :', data['accession'])
    print('Version    :', data['version'])
    print('Organismus :', data['organism'])
    print('Definition :', data['definition'])
    print('Länge      : {} nt'.format(len(seq)))
    print('GC-Gehalt  : {:.2f}%'.format(gc_content(seq)))
    print('Entropie   : {:.6f} bit/Base'.format(shannon_entropy(seq)))
    print('Basen      :', dict(Counter(seq)))

    print('\n--- Codon-Positionen (1,2,3) ---')
    for i, c in enumerate(base_counts_by_codon_position(seq), start=1):
        total = sum(c.values())
        freqs = {b: round(c.get(b, 0)/total, 4) for b in 'ACGT'} if total else {}
        print('Position {}: counts={} freqs={}'.format(i, dict(c), freqs))

    print('\n--- Häufigste 3-mer / Codons ---')
    print_top(Counter(codons(seq)), n=15, title='Top 15 Codons')

    print('\n--- Häufigste 4-mer ---')
    print_top(kmer_counts(seq, 4), n=15, title='Top 15 4-mer')

    print('\n--- CDS-Analyse ---')
    for feat in data['features']:
        analyze_feature(seq, feat)

    print('\n--- Übergangsmatrix A,C,G,T -> A,C,G,T ---')
    bases, M = transition_matrix(seq)
    print('Reihenfolge:', list(bases))
    print(M)

    print('\n--- Einfache algebraische Objekte ---')
    print('Vektor der Basenzahlen in Reihenfolge (A,C,G,T):')
    v = vector(ZZ, [seq.count('A'), seq.count('C'), seq.count('G'), seq.count('T')])
    print(v)
    print('Norm_1 =', sum(abs(x) for x in v))
    print('Norm_2^2 =', sum(x*x for x in v))

    print('\n--- Optional: de-Bruijn-artige Kanten für k=3 -> 4 ---')
    km3 = kmer_counts(seq, 3)
    edges = []
    for k4, cnt in kmer_counts(seq, 4).items():
        edges.append((k4[:3], k4[1:], cnt))
    for e in sorted(edges, key=lambda x: (-x[2], x[0], x[1]))[:20]:
        print('  {} -> {} : {}'.format(*e))


if __name__ == '__main__':
    main()

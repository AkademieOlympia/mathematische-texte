#!/usr/bin/env sage
# -*- coding: utf-8 -*-

"""
BM-nahe SageMath-Analyse einer GenBank-Datei (DNA.gb)
----------------------------------------------------
Diese Variante ergänzt die Standardanalyse um:
- modulo-12-Profil der Basenpositionen
- quaternionisch inspirierte Kodierung der Basen
- Fenster-Spektralanalyse für 3er-, 4er- und 12er-Raster
- Vergleich der Signaturen von ORF6 und ORF7a
- einfache Distanz- und Korrelationsmaße
- TSV-Exporte und PNG-Grafiken

Hinweis:
Die 'quaternionische' Kodierung ist hier eine diskrete Modellierung,
nicht eine biologisch etablierte Standarddarstellung.
"""

import os
import re
import math
from collections import Counter, defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sage.all import *

INPUT_FILE = 'DNA.gb'
OUTDIR = 'dna_bm_output'
os.makedirs(OUTDIR, exist_ok=True)

# -----------------------------
# Hilfsfunktionen
# -----------------------------

def read_text(path):
    with open(path, 'r', encoding='utf-8') as f:
        return f.read()


def extract_origin_sequence(text):
    m = re.search(r'ORIGIN(.*)//', text, re.S)
    if not m:
        raise ValueError('Kein ORIGIN-Block gefunden.')
    block = m.group(1)
    seq = ''.join(re.findall(r'[acgtunACGTUN]', block)).upper().replace('U', 'T')
    return seq


def parse_metadata(text):
    md = {}
    loc = re.search(r'^LOCUS\s+(\S+)\s+(\d+)\s+bp', text, re.M)
    if loc:
        md['locus'] = loc.group(1)
        md['length'] = int(loc.group(2))
    acc = re.search(r'^ACCESSION\s+(\S+)', text, re.M)
    if acc:
        md['accession'] = acc.group(1)
    org = re.search(r'\s+/organism="([^"]+)"', text)
    if org:
        md['organism'] = org.group(1)
    host = re.search(r'\s+/host="([^"]+)"', text)
    if host:
        md['host'] = host.group(1)
    country = re.search(r'\s+/country="([^"]+)"', text)
    if country:
        md['country'] = country.group(1)
    date = re.search(r'\s+/collection_date="([^"]+)"', text)
    if date:
        md['collection_date'] = date.group(1)
    iso = re.search(r'\s+/isolate="([^"]+)"', text)
    if iso:
        md['isolate'] = iso.group(1)
    source = re.search(r'\s+/isolation_source="([^"]+)"', text)
    if source:
        md['isolation_source'] = source.group(1)
    return md


def parse_cds_features(text):
    lines = text.splitlines()
    cds = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if re.match(r'^\s{5}CDS\s+', line):
            loc = line[21:].strip()
            m = re.match(r'(\d+)\.\.(\d+)', loc)
            if not m:
                i += 1
                continue
            start, end = int(m.group(1)), int(m.group(2))
            entry = {'start': start, 'end': end}
            i += 1
            while i < len(lines) and lines[i].startswith(' ' * 21 + '/'):
                qline = lines[i].strip()
                if qline.startswith('/gene='):
                    entry['gene'] = qline.split('=', 1)[1].strip().strip('"')
                elif qline.startswith('/product='):
                    entry['product'] = qline.split('=', 1)[1].strip().strip('"')
                elif qline.startswith('/translation='):
                    value = qline.split('=', 1)[1].strip()
                    if value.startswith('"') and value.endswith('"'):
                        entry['translation'] = value.strip('"')
                    else:
                        parts = [value.lstrip('"')]
                        i += 1
                        while i < len(lines):
                            t = lines[i].strip()
                            parts.append(t)
                            if t.endswith('"'):
                                break
                            i += 1
                        entry['translation'] = ''.join(parts).strip('"').replace(' ', '')
                i += 1
            cds.append(entry)
            continue
        i += 1
    return cds


def base_counts(seq):
    c = Counter(seq)
    return {b: c.get(b, 0) for b in 'ATCG'}


def gc_content(seq):
    c = base_counts(seq)
    return (c['G'] + c['C']) / len(seq) if seq else 0.0


def shannon_entropy(seq):
    n = len(seq)
    c = Counter(seq)
    H = 0.0
    for base in 'ATCG':
        p = c.get(base, 0) / n if n else 0
        if p > 0:
            H -= p * math.log(p, 2)
    return H


def codons(seq):
    return [seq[i:i+3] for i in range(0, len(seq) - len(seq) % 3, 3)]


def kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def translate(seq):
    aa = []
    for c in codons(seq):
        aa.append(CODON_TABLE.get(c, 'X'))
    return ''.join(aa)


def write_tsv(path, header, rows):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\t'.join(map(str, header)) + '\n')
        for row in rows:
            f.write('\t'.join(map(str, row)) + '\n')


# -----------------------------
# BM-nahe Kodierung
# -----------------------------
# Quaternionisch inspirierte Basiszuordnung:
# A ->  1
# T -> -1
# C ->  i
# G -> -i
# Dadurch bleibt die Abbildung 2-dimensional im Komplexen,
# aber algebraisch sauber für Fourier/Spektralanalyse nutzbar.

QMAP = {
    'A': CC(1),
    'T': CC(-1),
    'C': CC(I),
    'G': CC(-I),
}


def q_encode(seq):
    return vector(CC, [QMAP[b] for b in seq if b in QMAP])


def modulo_profile(seq, mod=12):
    prof = {b: [0]*mod for b in 'ATCG'}
    for idx, ch in enumerate(seq, start=1):
        if ch in prof:
            prof[ch][idx % mod] += 1
    return prof


def complex_profile_by_mod(seq, mod=12):
    vals = [CC(0)] * mod
    for idx, ch in enumerate(seq, start=1):
        if ch in QMAP:
            vals[idx % mod] += QMAP[ch]
    return vals


def dft_abs(values):
    # Diskrete Fourier-Analyse mit Sage/Komplexzahlen
    n = len(values)
    out = []
    for k in range(n):
        s = CC(0)
        for t, x in enumerate(values):
            s += x * exp(-2*pi*I*k*t/n)
        out.append(abs(complex(s)))
    return out


def window_sums(seq, w):
    enc = [QMAP[b] for b in seq if b in QMAP]
    out = []
    for i in range(0, len(enc) - w + 1):
        out.append(sum(enc[i:i+w]))
    return out


def spectral_signature(seq, w):
    ws = window_sums(seq, w)
    mags = [abs(complex(z)) for z in ws]
    if not mags:
        return {
            'window': w,
            'mean_abs': 0.0,
            'var_abs': 0.0,
            'max_abs': 0.0,
            'energy': 0.0,
        }
    mean_abs = sum(mags) / len(mags)
    var_abs = sum((x - mean_abs)**2 for x in mags) / len(mags)
    energy = sum(x*x for x in mags)
    return {
        'window': w,
        'mean_abs': mean_abs,
        'var_abs': var_abs,
        'max_abs': max(mags),
        'energy': energy,
    }


def segment_signature(seq):
    counts = base_counts(seq)
    gc = gc_content(seq)
    H = shannon_entropy(seq)
    qv = q_encode(seq)
    qsum = sum(qv) if len(qv) else CC(0)
    return {
        'length': len(seq),
        'A': counts['A'], 'T': counts['T'], 'C': counts['C'], 'G': counts['G'],
        'GC': gc,
        'entropy': H,
        'qsum_real': qsum.real(),
        'qsum_imag': qsum.imag(),
        'qsum_abs': abs(complex(qsum)),
    }


def cosine_similarity_real(v1, v2):
    # Für reelle Listen
    a = sum(x*x for x in v1) ** 0.5
    b = sum(x*x for x in v2) ** 0.5
    if a == 0 or b == 0:
        return 0.0
    return sum(x*y for x, y in zip(v1, v2)) / (a*b)


def flatten_mod_profile(prof):
    out = []
    for b in 'ATCG':
        out.extend(prof[b])
    return out


def plot_bar(values, labels, title, path):
    plt.figure(figsize=(10, 4))
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(values)), labels, rotation=90)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


def plot_line(values, title, path, xlabel='Index', ylabel='Wert'):
    plt.figure(figsize=(10, 4))
    plt.plot(range(len(values)), values)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


def plot_mod_profile(prof, title, path):
    plt.figure(figsize=(10, 5))
    xs = list(range(len(next(iter(prof.values())))))
    for b in 'ATCG':
        plt.plot(xs, prof[b], label=b)
    plt.title(title)
    plt.xlabel('Position mod 12')
    plt.ylabel('Anzahl')
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


# -----------------------------
# Hauptteil
# -----------------------------

def main():
    text = read_text(INPUT_FILE)
    seq = extract_origin_sequence(text)
    md = parse_metadata(text)
    cds_list = parse_cds_features(text)

    # Gesamtstatistik
    counts = base_counts(seq)
    gc = gc_content(seq)
    H = shannon_entropy(seq)
    all_codons = codons(seq)
    codon_count = Counter(all_codons)
    mer4_count = Counter(kmers(seq, 4))

    # Modulo-12 und komplexe Aggregation
    mod12 = modulo_profile(seq, 12)
    cmod12 = complex_profile_by_mod(seq, 12)
    cmod12_abs = [abs(complex(z)) for z in cmod12]
    cmod12_dft = dft_abs(cmod12)

    # globale Fenster-Signaturen
    sig3 = spectral_signature(seq, 3)
    sig4 = spectral_signature(seq, 4)
    sig12 = spectral_signature(seq, 12)

    # CDS-Segmente
    cds_rows = []
    cds_compare_rows = []
    cds_profiles = {}

    for entry in cds_list:
        gene = entry.get('gene', entry.get('product', 'CDS'))
        start, end = entry['start'], entry['end']
        sub = seq[start-1:end]
        aa = translate(sub)
        ann = entry.get('translation', '')
        ann_clean = ann.replace(' ', '').replace('\n', '')
        is_match = (aa == ann_clean) if ann_clean else ''

        s = segment_signature(sub)
        s3 = spectral_signature(sub, 3)
        s4 = spectral_signature(sub, 4)
        s12 = spectral_signature(sub, 12)
        prof = modulo_profile(sub, 12)
        cds_profiles[gene] = prof

        cds_rows.append([
            gene, start, end, len(sub),
            s['A'], s['T'], s['C'], s['G'],
            round(s['GC'], 6), round(s['entropy'], 6),
            round(float(s['qsum_real']), 6), round(float(s['qsum_imag']), 6), round(float(s['qsum_abs']), 6),
            round(s3['mean_abs'], 6), round(s3['energy'], 6),
            round(s4['mean_abs'], 6), round(s4['energy'], 6),
            round(s12['mean_abs'], 6), round(s12['energy'], 6),
            is_match
        ])

    # Vergleich ORF6 vs ORF7a, falls vorhanden
    if 'ORF6' in cds_profiles and 'ORF7a' in cds_profiles:
        f6 = flatten_mod_profile(cds_profiles['ORF6'])
        f7 = flatten_mod_profile(cds_profiles['ORF7a'])
        # Länge angleichen für kosinusbasierten Vergleich via Normierung getrennt je Vektor?
        # Hier gleiche Struktur (je 4*12 Einträge), also direkt vergleichbar.
        cos = cosine_similarity_real(f6, f7)
        dist = math.sqrt(sum((a-b)**2 for a, b in zip(f6, f7)))
        cds_compare_rows.append(['ORF6_vs_ORF7a', round(cos, 6), round(dist, 6)])

    # Summary schreiben
    with open(os.path.join(OUTDIR, 'summary.txt'), 'w', encoding='utf-8') as f:
        f.write('BM-nahe SageMath-Analyse der Datei DNA.gb\n')
        f.write('=====================================\n\n')
        for k, v in md.items():
            f.write(f'{k}: {v}\n')
        f.write(f'\nSequenzlänge: {len(seq)}\n')
        f.write(f'Basen: A={counts["A"]}, T={counts["T"]}, C={counts["C"]}, G={counts["G"]}\n')
        f.write(f'GC-Gehalt: {gc:.6f}\n')
        f.write(f'Shannon-Entropie: {H:.6f}\n')
        f.write('\nGlobale BM-nahe Signaturen:\n')
        for s in (sig3, sig4, sig12):
            f.write(f"Fenster {s['window']}: mean_abs={s['mean_abs']:.6f}, var_abs={s['var_abs']:.6f}, max_abs={s['max_abs']:.6f}, energy={s['energy']:.6f}\n")

    write_tsv(
        os.path.join(OUTDIR, 'base_counts.tsv'),
        ['Base', 'Count'],
        [[b, counts[b]] for b in 'ATCG']
    )

    write_tsv(
        os.path.join(OUTDIR, 'top_codons.tsv'),
        ['Codon', 'Count'],
        codon_count.most_common(32)
    )

    write_tsv(
        os.path.join(OUTDIR, 'top_4mers.tsv'),
        ['4mer', 'Count'],
        mer4_count.most_common(32)
    )

    # Mod12-Profile TSV
    rows = []
    for b in 'ATCG':
        for idx, val in enumerate(mod12[b]):
            rows.append([b, idx, val])
    write_tsv(os.path.join(OUTDIR, 'mod12_profile.tsv'), ['Base', 'ResidueClass', 'Count'], rows)

    write_tsv(
        os.path.join(OUTDIR, 'complex_mod12.tsv'),
        ['ResidueClass', 'Real', 'Imag', 'Abs'],
        [[i, round(float(z.real()), 6), round(float(z.imag()), 6), round(abs(complex(z)), 6)] for i, z in enumerate(cmod12)]
    )

    write_tsv(
        os.path.join(OUTDIR, 'complex_mod12_dft.tsv'),
        ['FrequencyIndex', 'Magnitude'],
        [[i, round(v, 6)] for i, v in enumerate(cmod12_dft)]
    )

    write_tsv(
        os.path.join(OUTDIR, 'cds_bm_signatures.tsv'),
        ['Gene', 'Start', 'End', 'Length', 'A', 'T', 'C', 'G', 'GC', 'Entropy',
         'QsumReal', 'QsumImag', 'QsumAbs',
         'MeanAbs_w3', 'Energy_w3', 'MeanAbs_w4', 'Energy_w4', 'MeanAbs_w12', 'Energy_w12', 'TranslationMatchesAnnotation'],
        cds_rows
    )

    if cds_compare_rows:
        write_tsv(
            os.path.join(OUTDIR, 'cds_comparison.tsv'),
            ['Comparison', 'CosineSimilarity_Mod12Profiles', 'EuclideanDistance_Mod12Profiles'],
            cds_compare_rows
        )

    # Grafiken
    top_cod = codon_count.most_common(20)
    plot_bar([v for _, v in top_cod], [k for k, _ in top_cod], 'Top-20 Codons', os.path.join(OUTDIR, 'top_codons.png'))

    top_m4 = mer4_count.most_common(20)
    plot_bar([v for _, v in top_m4], [k for k, _ in top_m4], 'Top-20 4-mere', os.path.join(OUTDIR, 'top_4mers.png'))

    plot_mod_profile(mod12, 'Modulo-12-Profil der Gesamtsequenz', os.path.join(OUTDIR, 'mod12_profile.png'))
    plot_line(cmod12_abs, 'Komplexe mod-12-Summen |z_r|', os.path.join(OUTDIR, 'complex_mod12_abs.png'), xlabel='Residuenklasse', ylabel='|z|')
    plot_line(cmod12_dft, 'DFT-Magnituden der komplexen mod-12-Summen', os.path.join(OUTDIR, 'complex_mod12_dft.png'), xlabel='Frequenzindex', ylabel='Magnitude')

    # CDS-spezifische Mod12-Plots
    for gene, prof in cds_profiles.items():
        plot_mod_profile(prof, f'Modulo-12-Profil {gene}', os.path.join(OUTDIR, f'mod12_{gene}.png'))

    print('Analyse abgeschlossen. Ergebnisse in:', OUTDIR)


if __name__ == '__main__':
    main()

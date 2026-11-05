#!/usr/bin/env python
"""Build the BioEmu-2 run set: FASTA sequences and system YAMLs.

Every mutant comes from the manuscript figures, not from a fresh scan, so that the
BioEmu-2 numbers sit directly beside the published CAMPARI and CALVADOS ones. Each
mutation is checked against the wild type before it is applied: if position 23 is not a
tryptophan, W23E is refused rather than silently written.

    python make_runs.py
"""
from __future__ import print_function

import os

HERE = os.path.dirname(os.path.abspath(__file__))
SEQ_DIR = os.path.join(HERE, 'sequences')
SYS_DIR = os.path.join(HERE, 'systems')

# ---------------------------------------------------------------- wild types

# p53 N-terminal transactivation domain, residues 1-61. Recovered from
# tests/fixtures/A39G_K24Y.fasta by reverting its two mutations.
P53_WT = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPD'

# PUMA BH3, from Figure 4B.
PUMA_WT = 'VEEEWAREIGAQLRRIADDLNAQYERRRQEEQHE'
PUMA_SCRAMBLE = 'ALAIEAWEQERQRIREDREHQYAVRGDQRLNEEE'

# UPF2 Q9HAU5, the 51-residue MIF4G-2 N-terminal disordered linker, and
# TDG Q13569 residues 1-50. Paste the sequences in to enable those two sets.
UPF2_WT = ''
TDG_WT = ''

# ---------------------------------------------------------------- mutant sets

# Figure 3F/G, the 20 designed double mutants. Strong attractive interactions found by
# the Interaction Map, disrupted by mutation to glutamic acid.
P53_DESIGNED = [
    'L35E_F54E', 'L26E_F54E', 'V31E_L45E', 'W23E_W53E', 'L25E_W53E',
    'L32E_G59E', 'W23E_F54E', 'P36E_F54E', 'P36E_W53E', 'L25E_F54E',
    'L32E_P58E', 'V31E_M44E', 'A39E_P58E', 'A39E_W53E', 'A39E_G59E',
    'L26E_W53E', 'A39E_F54E', 'P36E_P60E', 'A39E_P60E', 'M40E_G59E',
]

# Figure 3C, the four most frequent disease-associated single-point mutations.
P53_DISEASE = ['D49N', 'P4L', 'V31I', 'V10I']

# Figure 5F. UPF2 is repulsive, so these compact the ensemble.
UPF2_DESIGNED = [
    'E41Y_L42Y_E43Y', 'P35Y_D36Y_D37Y', 'L40Y_E46Y', 'V32I_E43I', 'P35Y_E43Y',
    'P35I_E41I', 'F3I', 'D36Y_E41Y', 'L45Y', 'V32K_L38E',
]

# Figure 6F. TDG is attractive, so these expand it.
TDG_DESIGNED = [
    'Q15E_A33E', 'S10E_T19E', 'S10E_P21E', 'S10E_Y18E', 'F17E_V35E',
    'G7E_T19E', 'Q15E_F20E', 'A16E_M32E', 'M26E_V35E', 'S10E_F17E',
]

GS_REPEATS = [4, 8, 16, 24, 32, 48, 64]


def apply_mutations(wt, label):
    """Apply an underscore-separated mutation label, checking each wild-type residue."""
    seq = list(wt)
    for token in label.split('_'):
        src, pos, dst = token[0], int(token[1:-1]), token[-1]
        if pos < 1 or pos > len(seq):
            raise ValueError('%s: position %d outside 1..%d' % (token, pos, len(seq)))
        if seq[pos - 1] != src:
            raise ValueError('%s: wild type has %s at %d, not %s'
                             % (token, seq[pos - 1], pos, src))
        seq[pos - 1] = dst
    return ''.join(seq)


def write_fasta(path, records):
    with open(path, 'w') as handle:
        for name, seq in records:
            handle.write('>%s\n%s\n' % (name, seq))
    return len(records)


def write_system(name, seq):
    with open(os.path.join(SYS_DIR, name + '.yaml'), 'w') as handle:
        handle.write('entities:\n')
        handle.write('  - entity_type: protein\n')
        handle.write('    asym_ids: [A]\n')
        handle.write('    sequence: %s\n' % seq)


def main():
    for d in (SEQ_DIR, SYS_DIR):
        if not os.path.isdir(d):
            os.makedirs(d)

    all_records = []

    # GS reference linkers
    gs = [('GS%d' % n, 'GS' * n) for n in GS_REPEATS]
    write_fasta(os.path.join(SEQ_DIR, 'gs_linkers.fasta'), gs)
    all_records += gs
    print('GS linkers      %d  (%s residues)'
          % (len(gs), ', '.join(str(len(s)) for _, s in gs)))

    # p53 wild type, designed mutants, disease mutants
    p53 = [('p53_WT', P53_WT)]
    for label in P53_DESIGNED:
        p53.append(('p53_' + label, apply_mutations(P53_WT, label)))
    write_fasta(os.path.join(SEQ_DIR, 'p53_designed.fasta'), p53)
    all_records += p53
    print('p53 designed    %d  (WT + %d mutants, all positions verified)'
          % (len(p53), len(P53_DESIGNED)))

    disease = [('p53_' + label, apply_mutations(P53_WT, label)) for label in P53_DISEASE]
    write_fasta(os.path.join(SEQ_DIR, 'p53_disease.fasta'), disease)
    all_records += disease
    print('p53 disease     %d  (positions verified)' % len(disease))

    # PUMA
    puma = [('PUMA_WT', PUMA_WT), ('PUMA_scramble', PUMA_SCRAMBLE)]
    write_fasta(os.path.join(SEQ_DIR, 'puma.fasta'), puma)
    all_records += puma
    same = sorted(PUMA_WT) == sorted(PUMA_SCRAMBLE)
    print('PUMA            %d  (WT %d aa, scramble %d aa, composition matches: %s)'
          % (len(puma), len(PUMA_WT), len(PUMA_SCRAMBLE), same))

    # UPF2 and TDG, only when their wild types have been pasted in
    for tag, wt, labels in (('upf2', UPF2_WT, UPF2_DESIGNED),
                            ('tdg', TDG_WT, TDG_DESIGNED)):
        if not wt:
            print('%-15s skipped  (paste the wild type into make_runs.py)' % tag)
            continue
        recs = [(tag + '_WT', wt)]
        for label in labels:
            recs.append((tag + '_' + label, apply_mutations(wt, label)))
        write_fasta(os.path.join(SEQ_DIR, tag + '.fasta'), recs)
        all_records += recs
        print('%-15s %d  (WT + %d mutants, positions verified)'
              % (tag, len(recs), len(labels)))

    for name, seq in all_records:
        write_system(name, seq)
    print()
    print('%d system yaml files in %s' % (len(all_records), SYS_DIR))


if __name__ == '__main__':
    main()

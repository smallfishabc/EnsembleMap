# BioEmu-2 run set

Sequences and system files for running the manuscript's designed mutants through BioEmu-2.

Every mutant here is taken from the paper's figures rather than from a fresh scan, so the
BioEmu-2 results land directly beside the published CAMPARI and CALVADOS numbers. Each
mutation was checked against its wild type before being written — all 24 p53 positions
carried the expected residue.

```
make_runs.py        builds sequences/ and systems/
sequences/*.fasta   the sequences, grouped
systems/*.yaml      one BioEmu-2 --system file per sequence
run_all.sh          the sampling loop
```

Regenerate with `python make_runs.py`.

## What this set is for

Figure 3F/G is the claim: of the 21 p53 sequences, **19 of 21 expanded** and **18 of 21
showed a decrease in attractive interactions** under CAMPARI. That is the design working.

The BioEmu-2 question is whether it reproduces those two counts. If it does, BioEmu-2 can
stand in for all-atom simulation in the design pipeline, and the screening step stops being
the expensive one. The designed mutants are the demonstration, which is why the bulk of the
runs are mutants rather than wild types — WT is here only as the baseline those two counts
are measured against, exactly as it appears as a point in Figure 3F.

## The runs

| group | n | what it is |
|---|---|---|
| GS linkers | 7 | GS4, GS8, GS16, GS24, GS32, GS48, GS64 — 8 to 128 residues |
| p53 designed | 21 | WT plus the 20 double mutants to glutamic acid from Fig. 3F/G |
| p53 disease | 4 | D49N, P4L, V31I, V10I, the four from Fig. 3C |
| PUMA | 2 | WT and the Fig. 4B scramble |
| UPF2 | 11 | Fig. 5F — needs the wild type pasted into `make_runs.py` |
| TDG | 11 | Fig. 6F — needs the wild type pasted into `make_runs.py` |

34 systems are built now; 56 once UPF2 and TDG have their wild types.

**GS first.** The map is `ln(p_1/p_0)` and `p_0` has to be fitted on BioEmu-2's own
linkers. Nothing else can be computed until that is done. The paper fits 8 to 32 repeats;
GS4 and GS48/GS64 extend the range at both ends, which also tests whether the exponent is
length-dependent in BioEmu-2 the way it is not in CAMPARI.

**p53 designed is the headline.** Twenty mutants, two positions each, hydrophobic or
aromatic residues mutated to glutamic acid to turn an attraction into a like-charge
repulsion.

**UPF2 and TDG go the other way**, which is why they are worth including. UPF2 is repulsive
and its mutants compact it (to Y and I); TDG is attractive and its mutants expand it (to E).
A model that only reproduces one direction has not reproduced the method.

## Two things to know

`UPF2_WT` and `TDG_WT` are empty in `make_runs.py`. UPF2 is the 51-residue MIF4G-2
N-terminal linker of Q9HAU5 and TDG is residues 1–50 of Q13569. Paste them in and rerun;
the mutation positions will be checked against them, so a wrong range will be caught rather
than silently mutated.

Both PUMA sequences are 34 residues and their compositions match exactly, which is
what a scramble should be. My first read of Figure 4B had four E residues in the
opening run rather than three; the composition check against the scramble located it.
`make_runs.py` reports the match on every run.

## Running

```bash
./run_all.sh            # all systems
./run_all.sh p53_       # just the p53 ones
```

Unsteered, fixed seed, one output directory per run. The comparison is between sequences,
so sampling noise should be shared.

## Analysis

BioEmu-2 is all-atom, so the published all-atom geometry applies: `scheme='ca'`, 0.8 nm,
minimum separation 3. Build the profile at run time rather than editing `profiles.py`:

```python
from dataclasses import replace
from ensemble_map import profiles
BIOEMU2 = replace(profiles.ALL_ATOM, name='bioemu2',
                  residue_index_offset=1,        # no ACE cap, unlike CAMPARI
                  reference_a=A, reference_b=B)  # from the GS fit
```

`residue_index_offset=1` matters — `ALL_ATOM` uses 0 because the CAMPARI topologies carry
an ACE cap at index 0. `check_residue_numbering(contact_df, sequence)` catches it.

Then per mutant, the two numbers Figure 3F/G reports: normalised ΔRg, and the normalised
change in the number of attractive interactions.

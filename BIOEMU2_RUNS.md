# BioEmu-2 run list

What to simulate with vanilla BioEmu-2, and which constrained runs are worth doing.
Protocols for how to analyse the output are in `BIOEMU2_PROTOCOL.md`.

Sequences are given as UniProt accession plus residue range rather than typed out, except
the GS linkers which are constructible exactly. The accessions are from memory and worth a
ten-second check before you pull them — one wrong residue invalidates a run.

---

## Part 1 — vanilla runs

All unsteered, so the ensembles are BioEmu-2's own:

```
bioemu2 sample-equilibrium --system <sys>.yaml --out runs/<name> \
    --num-samples <n> --batch-size 64 --sampler edm --seed 0
```

`examples/search_batch_size.py` sets the batch size for the card you have.

### A. Reference curve — do this first, nothing else works without it

The map is `ln(p_obs / p_ref)` and `p_ref` has to come from BioEmu-2's own linkers, not
from ABSINTH's. Four lengths so the exponent can be checked for length dependence.

| name | sequence | residues | samples |
|---|---|---|---|
| GS8 | `GSGSGSGSGSGSGSGS` | 16 | 2000 |
| GS16 | `GSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGS` | 32 | 2000 |
| GS24 | `(GS)×24` | 48 | 2000 |
| GS32 | `(GS)×32` | 64 | 2000 |

GS32 alone spans separations 3–63, which covers the published fitting range. Fit each
length separately first and compare `b` across them.

### B. Ground truth — the one comparison that licenses everything else

| name | system | samples | why |
|---|---|---|---|
| p53_TAD | P04637 residues 1–61 | 10000 | You can reproduce the published ABSINTH map on this pixel-identically. It is the only sequence where a wrong BioEmu-2 map is detectable. |
| E1A | P03255 residues 1–64 | 2000 | Your CG system, so it cross-checks against CALVADOS rather than ABSINTH. |

10,000 rather than 28,000: the map is a pairwise average over frames, so it converges well
before the ABSINTH ensemble size. Worth checking convergence by splitting the 10,000 in
half and comparing the two maps.

### C. Sensitivity — this is what demonstrates the method matters

A map that only reflects amino acid composition is a worse version of a hydrophobicity
plot. These pairs are the argument that it reflects patterning.

| name | system | samples |
|---|---|---|
| p53_TAD_scramble_1/2/3 | p53 1–61, composition preserved, order shuffled, three seeds | 2000 each |
| A39G_K24Y | the mutant in `tests/fixtures/` | 2000 |
| A39G_K24Y_WT | its wild type | 2000 |

Scrambles are already part of your methodology — the paper runs PUMA scrambles — so this
is continuous with it rather than a new control. The mutant pair is two point mutations,
which is the smallest sequence change that should still move the map.

### D. Known answers — independent validation

| name | system | samples | what is known |
|---|---|---|---|
| alpha_synuclein | P37840 residues 1–140 | 2000 | The long-range N-terminal/C-terminal contact is the most measured intramolecular interaction in the IDP field, from PRE and smFRET. If the map shows it without being told, that is a figure. |
| hnRNPA1_LCD | P09651 LC domain, ~186–320 | 2000 | Aromatic stickers drive the contacts, and the positions are known. |

### E. Negative control

| name | system | samples | expected |
|---|---|---|---|
| prothymosin_alpha | P06454 | 2000 | Highly charged and genuinely expanded. The map should be close to empty. |

A near-empty map here matters: it shows the normalization is not manufacturing structure
out of a polymer. The GS linkers do this too, but they are also the reference, so a
separate sequence is the honest version of the control.

**Fourteen runs in total.** A, then B, then the rest in any order.

---

## Part 2 — constrained and steered runs

### 1. Interface anchoring, unconstrained against constrained

This is `examples/sample_protein_protein_constraints.py` applied to your system. It already
runs the same complex twice, once with contacts and once without, same seed, same sampler —
which is the controlled comparison you want.

p53 TAD (chain A) with MDM2 (chain B, Q00987 roughly 25–109). Anchor the hydrophobic triad
into the cleft: F19, W23, L26. Their recipe for the bound is heavy atoms within 5 Å, add
1 Å, round up to the nearest 0.5 Å.

Three maps come out of it: chain A bound, chain A free, and the A–B interface.

### 2. Mutant against wild type, partner anchored

Same two chains, same contacts, same seed; only chain A's sequence changes.

```
ΔΔ(i,j) = [I_bound − I_free]_mutant − [I_bound − I_free]_WT
```

Every term is intra-chain with the same reference curve, so the reference cancels. This is
the readout that says which intramolecular interactions a partner has to break, and which
mutations change that. Take the top candidates from `mutation_scan.design_mutants()` so the
mutants are the ones your method already proposes.

Ten mutants plus WT across two partners is 22 complex runs plus 11 free runs.

### 3. Map transfer — steer BioEmu-2 with the ABSINTH map

`U(x) = −Σ I_ij · 1[d_ij < c]` over the published p53 pairs, as a `Potential` subclass.
`examples/tilt_a_sample_distribution_with_steering.py` is the template, and its own
verification is the control you need: sample unsteered, reweight those frames by `exp(−U)`,
and compare the reweighted distribution against the steered one. If they agree, the sampler
is targeting `p_base·exp[−U]` and the steering is doing what it claims.

`FKCSampler(potentials=[...], ess_threshold=0.8, weight_update="guidance", reheat="sigma")`
— the `reheat="sigma"` is required with `"guidance"`.

### 4. Contact dose–response

The one I would most like to see. Constrain only the top 1, then 2, 3, 5, 10 attractive
pairs from the published map, and after each run compute the *whole* map and measure how
much of the rest appeared on its own.

That gives a curve: how much of an Ensemble Map is determined by its strongest few
contacts. If five pairs reproduce most of the map, the map has a small number of degrees of
freedom and those pairs are the ones worth mutating. If it stays flat, every arc is
carrying independent information. Either answer is a result, and it is five cheap runs.

### 5. Scramble with the same constraints

Apply the p53 contacts to a scrambled p53 sequence. If the map still forms, the constraints
are doing the work; if it does not, the sequence is. This is the control that keeps run 3
from being circular.

---

## Settings to hold fixed across everything

- `--sampler edm` for all vanilla runs, so nothing is biased
- the same `--seed` for every run in a comparison — the comparison is between sequences, so
  sampling noise should be shared
- query-only MSA where the example uses one, which avoids alignment infrastructure entirely
- one output directory per run; the CLI refuses to write into an existing one

## Analysis side

BioEmu-2 is all-atom, so use the all-atom geometry: `scheme='ca'`, cutoff 0.8 nm, min
separation 3. Build the profile at run time rather than editing `profiles.py`:

```python
from dataclasses import replace
from ensemble_map import profiles
BIOEMU2 = replace(profiles.ALL_ATOM, name='bioemu2',
                  residue_index_offset=1,        # BioEmu-2 builds from sequence, no caps
                  reference_a=A, reference_b=B)  # from part 1A
```

`residue_index_offset=1` matters: `ALL_ATOM` uses 0 because the CAMPARI topologies carry an
ACE cap at index 0. Check it with `check_residue_numbering(contact_df, sequence)` before
trusting any map.

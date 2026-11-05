# BioEmu-2 protocols

Three ways to use BioEmu-2 with the Ensemble Map. Each section assumes BioEmu-2 is already
installed and working; installation is set aside in the last section.

Nothing here requires a change to `src/ensemble_map/`. Where a new reference curve or a new
residue offset is needed, it is built at run time with `dataclasses.replace`.

- [Shared setup](#shared-setup)
- [Protocol 1 — calibrate against BioEmu-2](#protocol-1--calibrate-against-bioemu-2)
- [Protocol 2 — the map as a steering potential](#protocol-2--the-map-as-a-steering-potential)
- [Protocol 3 — mutant against wild type with a partner](#protocol-3--mutant-against-wild-type-with-a-partner)
- [Installation, set aside](#installation-set-aside)

---

## Shared setup

### The handoff

Two environments that never meet, exchanging files.

| | BioEmu-2 | Ensemble Map |
|---|---|---|
| Python | 3.12+ | 3.10 (`.conda-env`) |
| Platform | Linux + NVIDIA GPU | anywhere |
| Writes | `samples.xtc`, `topology.cif`, `physicality.json` | `*_interaction_*.csv`, figures |
| Reads | `.cif`, an `mdtraj.Trajectory` | `.cif`, `.xtc` |

mdtraj 1.9.9 registers `.cif` and `.pdbx`, so `topology.cif` loads directly and no converter
is needed. `propagate()` takes an `mdtraj.Trajectory`, so the path back is open too.

### Geometry

BioEmu-2 is all-atom, so it lands on the published all-atom geometry: Cα distances, 0.8 nm
cutoff, minimum separation 3. Use `--profile all_atom` for anything that only needs the
geometry.

### Residue numbering — the one thing to get right first

`ALL_ATOM` carries `residue_index_offset=0` because the CAMPARI topologies have an ACE cap
at residue index 0, which makes raw mdtraj indices already 1-based. **BioEmu-2 builds from
sequence and has no caps**, so its residue index 0 *is* sequence position 1 and the offset
must be 1, as for the coarse-grained profile.

Build the profile at run time rather than editing `profiles.py`:

```python
from dataclasses import replace
from ensemble_map import profiles

BIOEMU2 = replace(profiles.ALL_ATOM,
                  name='bioemu2',
                  residue_index_offset=1,       # no terminal caps
                  reference_a=A, reference_b=B) # from Protocol 1
```

`profiles.get_profile()` accepts a `Profile` object, so this passes anywhere a profile name
does. Confirm the numbering before trusting any map:

```python
from ensemble_map.contact_map_generation import check_residue_numbering
check_residue_numbering(contact_df, sequence)   # must return True, prints nothing
```

It returns `False` and prints the span it found when `r_1`/`r_2` do not run `1..len(sequence)`.

### Command spellings

The `ensemble-map` and `calibrate` options below are read from the source and are exact. The
BioEmu-2 flags follow its docs, but confirm the spellings once with
`bioemu2 sample-equilibrium --help` and `bioemu2 propagate --help` before scripting a sweep.

### Naming

One directory per run, since the BioEmu-2 CLI refuses to write into an existing output
directory:

```
runs/<system>_<variant>_<n-samples>_<seed>/
    system.yaml  samples.xtc  topology.cif  physicality.json
```

---

## Protocol 1 — calibrate against BioEmu-2

**Purpose.** Fit the reference curve `p_ref(s) = a·s^b` on BioEmu-2's own glycine-serine
linkers, so BioEmu-2 maps are measured against BioEmu-2's polymer baseline and are directly
comparable with the published ABSINTH maps. Also establishes how far BioEmu-2's chain
statistics sit from ABSINTH's.

### 1.1 Reference systems

Four GS linkers, `(GS)n` for n = 8, 16, 24, 32, giving 16, 32, 48 and 64 residues. The
longest spans separations 3–63, which covers the published fitting range.

```yaml
# runs/gs32/system.yaml
entities:
  - entity_type: protein
    asym_ids: [A]
    sequence: GSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGS
```

### 1.2 Sample

```bash
bioemu2 sample-equilibrium \
    --system runs/gs32/system.yaml \
    --out    runs/gs32 \
    --num-samples 2000 \
    --batch-size 64 \
    --sampler edm \
    --seed 0
```

Unsteered (`edm`) is correct here — the reference must be the model's own unbiased polymer
behaviour. 2000 conformations is enough for a separation-averaged curve, since every chain
contributes many pairs at each separation.

### 1.3 Fit

```bash
python -m ensemble_map.util.calibrate \
    runs/gs32/topology.cif runs/gs32/samples.xtc \
    --cutoff 0.8 --sep 3 60 \
    --source "bioemu2 (GS)32, 2000 samples, edm, seed 0" \
    --out references/bioemu2_gs32.json
```

Repeat for each length. The fit records `r_squared`, `residuals`, `a_stderr` and `b_stderr`,
so read them rather than just `a` and `b`.

**Checks, in order of what they tell you:**

1. `r_squared` on each length — does a power law describe BioEmu-2's linkers at all.
2. `b` across the four lengths — should agree within a couple of standard errors. If short
   and long linkers disagree, BioEmu-2's scaling is length-dependent and the pooled fit is
   the wrong summary.
3. `b` against the published `-1.33` — how far BioEmu-2's chain statistics are from ABSINTH.

For a pooled fit across lengths, `calibrate` takes one topology, so pool at the library
level instead:

```python
from ensemble_map.util.calibrate import (pair_separations_and_distances,
                                         contact_probability_by_separation, fit_power_law)
# accumulate per-separation probabilities from each length, average, then fit once
```

### 1.4 Target ensemble

Same thing for the real sequence, keeping the sampler unsteered:

```bash
bioemu2 sample-equilibrium --system runs/p53_wt/system.yaml \
    --out runs/p53_wt --num-samples 28000 --batch-size 64 --sampler edm --seed 0
```

28,000 matches the published ABSINTH ensemble size, which makes the comparison like for
like. Run it in several seeded batches if the session limit gets in the way, then pass all
the `samples.xtc` files to one analysis.

### 1.5 Map and compare

```bash
ensemble-map --pdb runs/p53_wt/topology.cif \
             --xtc runs/p53_wt/samples.xtc \
             --name p53_bioemu2 --profile all_atom --single_traj --no_show \
             --out runs/p53_wt/map
```

That uses the published `a, b`. To put it on BioEmu-2's own reference, re-normalize with the
`BIOEMU2` profile from the shared section:

```python
from ensemble_map import normalization
mapped = normalization.normalization(contact_df, profile=BIOEMU2)
```

**Compare against the ABSINTH map** on three things, not on pixels — two force fields will
never be pixel-identical and that is not the question:

| measure | what it answers |
|---|---|
| Spearman ρ of `relative_strength` over shared pairs | does BioEmu-2 rank interactions the same way |
| 5×5 confusion matrix of `plot_value` | where the levels disagree, and in which direction |
| per-level counts | whether BioEmu-2 is systematically more or less contact-rich |

### 1.6 The calibration knob

MSA depth is the main lever on ensemble breadth: shallower alignments give broader, more
flexible ensembles. Sweep it and pick the depth whose map best matches ABSINTH:

```bash
for depth in 8 32 128 512 full; do
    bioemu2 sample-equilibrium --system runs/p53_wt/system.yaml \
        --out runs/p53_msa_${depth} --num-samples 2000 \
        --msa-subsample ${depth} --sampler edm --seed 0
done
```

Report the depth that matches. That single number is the calibration result, and it is
reusable for every later run on this system.

**Deliverable.** `references/bioemu2_*.json`, a `ρ` and a confusion matrix against the
published p53 map, and one chosen MSA depth.

---

## Protocol 2 — the map as a steering potential

**Purpose.** Push BioEmu-2's generation toward a measured interaction map, so the
information in an expensive ABSINTH ensemble is carried into a cheap generator.

### 2.1 Export the target

From the published interaction CSV, keep the pairs that carry signal:

```python
import pandas as pd
m = pd.read_csv('All_atom_p53_interaction/.../interaction_0317.csv')
target = m[m['plot_value'].abs() >= 1][['r_1', 'r_2', 'relative_strength', 'p_sem']]
```

`relative_strength` is `I = ln(p_obs/p_ref)`, which is the quantity that belongs in the
energy. `p_sem` is the per-pair standard error and makes a natural confidence weight, so
well-sampled pairs pull harder than marginal ones.

### 2.2 The potential

`Potential` wants `prepare(system)` and `__call__(coords_angstrom, *, t)` returning one
energy per structure. Coordinates arrive in ångström, so the 0.8 nm cutoff is 8 Å. Guidance
needs differentiability, hence the sigmoid in place of a hard contact indicator; distances
are rotation and translation invariant, which is what the sampler requires.

```python
import torch
from bioemu2.denoising.steering.potentials import Potential

class MapPotential(Potential):
    """Energy is low when the ensemble reproduces the target interaction map."""

    def __init__(self, pairs, strengths, cutoff_ang=8.0, tau_ang=1.0, parameters=None):
        super().__init__(parameters=parameters or {})
        self.pairs = pairs            # (N, 2), 1-based sequence positions
        self.I = torch.as_tensor(strengths, dtype=torch.float32)
        self.c, self.tau = cutoff_ang, tau_ang

    def prepare(self, system):
        # 1-based sequence position -> CA atom index, chain A
        self.i = torch.as_tensor([ca_index(system, 'A', r) for r in self.pairs[:, 0]])
        self.j = torch.as_tensor([ca_index(system, 'A', r) for r in self.pairs[:, 1]])

    def __call__(self, coords_angstrom, *, t):
        d = (coords_angstrom[..., self.i, :] - coords_angstrom[..., self.j, :]).norm(dim=-1)
        w = torch.sigmoid((self.c - d) / self.tau)      # soft contact indicator
        return -(self.I.to(d) * w).sum(-1)              # one energy per structure
```

Attractive pairs have `I > 0` and lower the energy when in contact; repulsive pairs have
`I < 0` and raise it. No separate sign handling is needed.

### 2.3 Run with guidance

```python
from bioemu2 import EmuSession, GuidedFKSampler, System, load_model
from bioemu2.denoising.steering.potentials import get_default_potentials

system = System.from_yaml('runs/p53_wt/system.yaml')
emu = EmuSession(system, model=load_model(), mode='equilibrium')

samples = emu.sample(
    num_samples=2000, batch_size=64,
    sampler=GuidedFKSampler(
        potentials=get_default_potentials() + [
            MapPotential(pairs, strengths, parameters={'guidance_weight': 0.05})],
        resampling_interval=None,     # guidance only to start
        num_gd_steps=20),
    seed=0)
```

Two things that will silently cost you a run if missed:

- **`guidance_weight` defaults to 0**, which means no guidance at all. It must be set.
- **Passing `potentials` replaces the defaults.** Concatenate `get_default_potentials()` or
  the physicality terms are gone and the structures degrade.

For a target distribution rather than a nudge, `FKCSampler(weight_update="smc")` targets
`p_base(x)·exp[−U(x)]` and needs only energy *values*, so the hard contact indicator can be
used directly and the sigmoid is unnecessary.

### 2.4 Weight sweep

```bash
for w in 0 0.01 0.05 0.2 1.0; do   # run the sweep through the Python API above
    ...
done
```

For each weight, recompute the map and record two numbers: agreement with the target, and
the physicality score from `physicality.json`. Agreement rises and physicality falls; take
the knee. Report the curve — it is the evidence that the steering is doing something real
rather than crushing the ensemble.

### 2.5 The control that makes this publishable

Steering toward a map and then recovering that map is circular. Hold pairs out:

- **Split-half.** Steer with a random half of the pairs, score agreement on the other half.
- **Sign-blind.** Steer with the attractive levels only (+2, +1) and ask whether the
  repulsive structure (−2) appears without being supplied.
- **Transfer.** Steer with the p53 map, then sample a *different* sequence with no potential
  and confirm the map does not persist — i.e. the potential acts on this system, not on the
  model globally.

The split-half number is the one to quote.

### 2.6 Cheaper first passes

Before writing any potential, two routes that need no new code at all:

1. **`ContactPotential`** on the level +2 pairs. It penalises pairs beyond a contact
   distance, which is the attractive half of the map.
2. **Contact conditioning** in the system YAML, which is a model input rather than an
   energy:

```yaml
contacts:
  - token1: {asym_id: A, res_id: 19}
    token2: {asym_id: A, res_id: 48}
    max_distance: 8.0
```

Conditioning is binary, so to express a graded map, draw a contact subset per batch —
include pair `(i,j)` with probability `cont_prob(i,j)` — run several batches with different
draws, and pool. The pooled ensemble reproduces the graded map in expectation.

**Deliverable.** An agreement-versus-physicality curve, a chosen `guidance_weight`, and a
split-half score.

---

## Protocol 3 — mutant against wild type with a partner

**Purpose.** Predict which mutations change partner binding, and which intramolecular
interactions a partner has to break. This is the one that extends the existing mutation scan
rather than re-validating it, and it is cheap because BioEmu-2 takes sequences — a mutant is
one changed letter, with no rebuilding and no re-equilibration.

### 3.1 Choose the mutants

From the existing map, using the shipped scan so the candidates are the ones the method
already proposes:

```python
from ensemble_map.util import mutation_scan
designs = mutation_scan.scan(interaction, sequence, name='p53',
                             interaction_type='attractive')
```

`select_pairs(interaction, contact_probability_cutoff=0.01, min_separation=5)` and
`design_mutants(...)` are the underlying steps if you want to filter differently. Take the
top 10.

### 3.2 Build the systems

Per sequence, two systems. Free:

```yaml
# runs/p53_wt_free/system.yaml
entities:
  - entity_type: protein
    asym_ids: [A]
    sequence: <p53 TAD 1-61>
```

Bound, with the interface anchored so sampling stays associated instead of drifting apart:

```yaml
# runs/p53_wt_mdm2/system.yaml
entities:
  - entity_type: protein
    asym_ids: [A]
    sequence: <p53 TAD 1-61>
  - entity_type: protein
    asym_ids: [B]
    sequence: <MDM2 25-109>
contacts:
  - token1: {asym_id: A, res_id: 23}    # p53 W23
    token2: {asym_id: B, res_id: 30}    # into the MDM2 cleft
    max_distance: 8.0
  - token1: {asym_id: A, res_id: 19}    # p53 F19
    token2: {asym_id: B, res_id: 37}
    max_distance: 8.0
```

`res_id` is 1-based **within the sequence as written in the YAML**, not UniProt numbering.
If chain A starts at p53 residue 1 then W23 is `res_id: 23`; if the construct is truncated,
subtract the offset. Fix the partner's construct boundaries once and keep them identical
across every mutant, or the comparison is not controlled.

### 3.3 Sample

```bash
for v in wt mut01 mut02 ... mut10; do
  for s in free mdm2 taz1 taz2; do
    bioemu2 sample-equilibrium \
        --system runs/p53_${v}_${s}/system.yaml \
        --out    runs/p53_${v}_${s} \
        --num-samples 2000 --batch-size 64 --sampler edm --seed 0
  done
done
```

11 variants × 4 systems = 44 runs at ~161 tokens each. Time one run first, then multiply;
`examples/search_batch_size.py` sets the batch size for the hardware you land on.

Keep the seed fixed across variants. The comparison is between sequences, so the sampling
noise should be shared.

### 3.4 Split the contact table by chain

The contact table already carries `r_1_chain` and `r_2_chain`, so the split needs no new
code:

```python
intra_A    = contact_df[(contact_df.r_1_chain == 0) & (contact_df.r_2_chain == 0)]
interface  = contact_df[contact_df.r_1_chain != contact_df.r_2_chain]
```

**`intra_A` is immediately valid.** Both residues are in chain A, so `distance = r_2 - r_1`
is a true sequence separation and `a·s^b` normalizes it correctly. Run the existing
`normalization()` on this subset and nothing else is required.

**`interface` needs its own reference**, because `a·s^b` is a function of sequence
separation and two chains do not have one. Use the marginal product of the contact table —
pointwise mutual information over the A×B block:

```python
import numpy as np
P = interface.pivot(index='r_1', columns='r_2', values='cont_prob').fillna(0).values
expected = P.sum(1, keepdims=True) * P.sum(0, keepdims=True) / P.sum()
I_inter  = np.log(np.where(P > 0, P, np.nan) / expected)
```

This does for two chains what the polymer curve does for one. The polymer curve removes
propensity that comes from sequence separation; the marginal product removes propensity that
comes from each residue's overall stickiness and each position's accessibility. What is left
is specific pairing: positive means this pair contacts more than either partner's general
tendency explains. It keeps the form `ln(p_obs/p_ref)`, so the five levels and the arc
drawing work unchanged, with only the thresholds wanting their own calibration.

Note that chain B's residues come out in continuous numbering — with chain A at 1–61, chain
B is 62 onward — so subtract 61 to report B in its own numbering.

As a physical check on the statistical null, sample chain A against several **scrambled**
partner sequences and build the reference indexed by residue-*type* pair, a 20×20 table,
since scrambling moves identities off their positions. Agreement with the marginal product
confirms the independence assumption; disagreement means the scramble is the one to trust.

### 3.5 The quantity

```
ΔΔ(i,j) = [ I_bound(i,j) − I_free(i,j) ]_mutant − [ I_bound(i,j) − I_free(i,j) ]_WT
```

Every term is intra-chain A, every term uses the same GS reference curve, and the reference
cancels in the difference. **So this needs none of the interface machinery above** — the
interface map is a separate deliverable, not a prerequisite. Protocol 3 can start as soon as
Protocol 1 gives a reference curve.

Per-mutant scalar, aggregated the same way the method already aggregates, so the numbers stay
comparable with the published scan:

```python
from ensemble_map.util.mutation_scan import chunk_strength, strength_lookup
lookup = strength_lookup(delta_delta)
score  = chunk_strength(lookup, r1, r2, window=1)    # Eq. 3 weighting
```

### 3.6 Read-outs

| comparison | what it shows |
|---|---|
| `I_bound − I_free`, WT | which intramolecular interactions the partner breaks |
| `ΔΔ` ranked over mutants | which mutations change binding, and in which direction |
| `ΔΔ` across partners | shared hub positions against partner-specific ones |
| `ΔΔ` rank against the scan's rank | whether the intramolecular scan predicts binding effects |

The last row is the result either way. Agreement says the existing scan already predicts
partner effects; disagreement locates the effects it cannot see, which is the argument for
doing the bound-state calculation at all.

### 3.7 Bound-state relaxation

For time-resolved behaviour within the bound state, propagate from the sampled complexes:

```bash
bioemu2 propagate --x0-cif runs/p53_wt_mdm2/topology.cif \
                  --x0-xtc runs/p53_wt_mdm2/samples.xtc \
                  --lag-time-ns 10 --out runs/p53_wt_mdm2_lag10
```

Each output is an endpoint at that lag, not an intermediate frame, so a series of lags gives
a series of ensembles; `simulate_trajectory()` chains steps when time-ordered frames are
wanted. Compute a map per lag and the ordering of contact formation is readable from the
sequence of maps. Resampling samplers cannot be used for multi-step trajectories, so keep
`edm` or `edm_gd` here.

### 3.8 Good first subject

E1A is the stronger multi-partner case than p53 — it binds CBP/p300, pRb, CtBP and TBP, and
`E1A_1_64` is already the system in `.vscode/launch.json`. Only chain B changes between
partners, so the marginal cost of a fourth partner is one more sequence.

**Deliverable.** A ranked `ΔΔ` table over mutants and partners, the WT `I_bound − I_free`
map showing what each partner breaks, and the rank correlation against the existing scan.

---

## Installation, set aside

Not needed to read the protocols above. Kept here so the decisions are in one place when we
come back to it.

### Where

No hosted option exists — the repository has no Azure guide, no Dockerfile and no managed
endpoint. Linux with an NVIDIA GPU, locally or rented.

| | GPU | Session | Cost |
|---|---|---|---|
| Google Colab | T4, 15 GB | 12 h, 90 min idle | free / Pro |
| Kaggle | T4 ×2 | 12 h, 30 h per week | free |
| WSL2 on this machine | whatever card is in it | unlimited | free |
| University cluster | A100 / H100 | batch | free |
| RunPod / Vast / Lambda | T4 to A100 | unlimited | ~$0.20–0.40/h T4, ~$1–2/h A100 |

Supported cards for the `cuda` extra are compute capability 7.5 and above — T4, A100,
RTX A6000, L40S, H100, B200. V100 (7.0) needs the `cuda12` extra instead. Below 7.0 is
unsupported, which rules out Kaggle's P100 option: select T4 there.

Colab's free T4 is exactly 7.5 with 15 GB, and at ~161 tokens none of the memory scaling
matters. Session limits are the real constraint for Protocol 3's 44 runs, which is when a
rented A100 earns its cost.

### How

The repository ships `uv.lock` and `.python-version`, so let `uv` supply Python 3.12 rather
than fighting whatever the host has:

```bash
nvidia-smi --query-gpu=name,compute_cap --format=csv   # check the card first
pip install uv
uv sync --extra cuda                                   # or --extra cuda12 on a V100
```

Three points of friction, none of them hard:

- **Two environments, permanently.** BioEmu-2 needs Python ≥3.12 and this package is pinned
  to 3.10 with mdtraj 1.9.9. They never share an environment; the handoff is files, which
  works in both directions.
- **MSAs** come from ColabFold retrieval or your own A3M files. Protocol 1.6 treats the
  depth as a parameter, so decide where they come from before that sweep.
- **Examples are `.py` with `# %%` cell markers**, not notebooks. Jupytext converts them, or
  the cells paste into `Google Colab/EnsembleMap_2026.ipynb`.

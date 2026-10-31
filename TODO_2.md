# TODO 2

Everything outstanding that is not in `TODO.md`. `TODO.md` holds the five things that
block releasing alongside the paper; this holds the rest, which is larger but none of it
urgent.

## Paper

- [ ] **Align the manuscript with the new name.** The title, the abstract, the Fig. 3A and
      3B captions and the Methods heading all say "Interaction Map". The code now says
      Ensemble Map. One or the other has to move, and the manuscript is the one still open.

- [ ] Decide whether the quantity keeps its name. Eq. 2 defines `I = ln(p_obs / p_ref)` as
      the **interaction strength**, and the code still calls it that: `interaction_plot.py`,
      `interaction_map()`, `relative_strength`, `plot_value`. That reads consistently — the
      product is the Ensemble Map, the thing it draws is an interaction — but it is worth
      being deliberate about rather than inheriting.

## Colab and Binder

- [ ] **Run `Google Colab/EnsembleMap_2026.ipynb` end to end on a GPU.** Only its analysis
      half has been run, against the committed fixture. The simulation half has never been
      executed in this project.

- [ ] Decide how the notebook is hosted. Opening it in Colab makes a copy, which then
      drifts. The usual arrangement is to keep the notebook in the repository and link to
      it with an "Open in Colab" badge, so there is one version.

- [ ] Update the notebook's install line once the branch is merged. It currently installs
      from `@refactor_20261024`; it should install from a tag.


- [ ] **Binder**, if it is wanted for reproducing the analysis without a GPU: add
      `jupyterlab` and `ipykernel` to `environment.yml`, add a `postBuild` file containing
      `pip install -e . --no-deps`, then add the badge. `tutorials/calibration.ipynb`
      already runs entirely on the committed fixture, so it works as the landing notebook.

## Analyses still in notebooks

Described in the paper, not yet in the package. Sources are under `temp_examples/`.

- [ ] **PCA clustering** of conformations by interaction pattern (Fig. 3D, E). Filters to
      separation 4–18 and attractive pairs, builds a per-frame binary contact matrix, then
      `PCA(n_components=2)` with `KMeans(n_clusters=5)`.

- [ ] **Solution-sensitivity regression** against attractive interaction count
      (Fig. 4D, E). Blocked on an external `simulationdatabase` package that is not in this
      repository; decide whether it comes along or the regression stays a notebook.

- [ ] Hydrophobic scaling per chunk (Fig. 5C/D, 6C/D) is **partly done** —
      `util/mutation_scan.py` computes per-chunk Kyte-Doolittle, charge and aromatic
      counts. What is missing is the histogram figure.

## Still on old branches

- [ ] Chunk aggregation and chunk plotting (`chunkreal`). Note this is a *different*
      quantity from the Eq. 3 window now in `util/mutation_scan.py`: it averages five
      residues with equal weight and has its own reference curve, a = 1.29, b = −1.22.
      Decide whether it is still part of the method.

- [ ] The mutation-target sweep (`clean_up_structure:pick_ensemble.py`), which ranks
      contacts by how much the ensemble changes when they are formed.

- [ ] `interaction_strength.py` cannot run on any branch — it references undefined names.
      Its only real use was counting interactions, which the regression above needs.
      Delete it or rewrite it as one function.

## Decisions to settle

- [ ] **Which contact-probability cutoff is the published one.** 0.1 with a minimum
      separation of 10 in the 2024 scan notebook; 0.01 with a minimum separation of 5 in
      the 2025 Colab notebook; none in the shipped analysis. They are parameters now, but
      the paper reports one pair. Measured effect below.

- [ ] **Which interaction thresholds are the published ones.** The current pairwise set
      reproduces the published p53 figure, so that is settled. Whether the chunk set is a
      separate calibrated choice is not. The four sets in circulation are below.

### Measured effect of the p₁ floor

On the published p53 map, 28,000 frames. The floor touches only the repulsive side —
the attractive counts do not move at all.

| p₁ floor | ≈ frames | lvl 2 | lvl 1 | lvl 0 | lvl −1 | lvl −2 |
|---|---|---|---|---|---|---|
| 0 (published) | 0 | **60** | **129** | 625 | 166 | **731** |
| 1/28000 | 1 | 60 | 129 | 696 | 166 | 660 |
| 0.001 | 28 | 60 | 129 | 1101 | 166 | **255** |
| 0.01 | 280 | 60 | 129 | 1406 | 85 | 31 |

So it is a counting decision rather than a drawing one. Between a floor of 0 and 0.001 the
number of repulsive arcs drawn falls from 897 to 421, but only 5% of pixels change and
0.16% change visibly, because the arcs removed are drawn at near-zero width — the width
scales with contact probability. Anything that tallies repulsive interactions is sensitive
to this; the figure is not.

### The four threshold sets

The same p53 data and the same reference curve under each set, to show the size of the
choice.

| set | lvl 2 | lvl 1 | lvl 0 | lvl −1 | lvl −2 |
|---|---|---|---|---|---|
| pairwise, current (0811) | 60 | 129 | 625 | 166 | 731 |
| pairwise, older (1115 / chunkreal) | 60 | 129 | 625 | 87 | 810 |
| chunk (chunkreal) | 115 | 74 | 547 | 78 | 897 |
| CG edges on all-atom data | 28 | 87 | 699 | 166 | 731 |

## Packaging, only if publishing to PyPI

- [ ] The dependency pins are exact `==`, which will downgrade a user's numpy and pandas on
      install, and mdtraj 1.9.9 has no wheels so pip will try to compile it. For a package
      people install, those want to be ranges, with `environment.yml` keeping the exact
      versions for reproducing the paper. This conflicts with staying on mdtraj 1.9.x, so
      it may be better not to publish to PyPI at all and to point people at
      `environment.yml` instead.

- [ ] If publishing: upload to TestPyPI first and install from it into a clean environment.
      A bad PyPI upload cannot be replaced, only yanked.

## Small

- [ ] `TUTORIAL.md` predates mutation scanning and calibration, so it covers the map and
      the conditional analysis only. Every command in it has been run and still passes.

- [ ] The package can be imported without mdtraj — the map, the reference curve and the
      figure need only numpy, pandas and matplotlib, and the trajectory readers load on
      first use. This was added for a reason that turned out not to apply. It is tested and
      harmless; revert it if the smaller change surface is worth more.

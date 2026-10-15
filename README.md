# gammaLoss

Tables, example configurations, analysis code and data for the `gammaLoss`
thermal photon module of C-SCAPE (thermal photon emission from the QGP and
hadron gas, and reabsorption of thermal and jet photons in the medium).

The module itself lives in C-SCAPE (`src/photons`). This repository holds what
it needs at run time and what is used to run and analyze photon simulations.

## Installation

From the C-SCAPE `external_packages` directory:

```bash
./get_gammaLoss.sh
```

This clones the repository into `external_packages/gammaLoss`. Then build
C-SCAPE with ROOT (MathMore is required for the photon modules):

```bash
cd build
cmake .. -DUSE_ROOT=ON -DUSE_MUSIC=ON -DUSE_ISS=ON -DUSE_SMASH=ON
make -j
```

With the repository present, the analysis programs in `analysis/` are built
as well.

## Contents

| directory | contents |
|-----------|----------|
| `tables/` | Rapp HRG photon rate and pion/kaon chemical potential tables used by the hadron gas rate (`<hadronGasTablePath>`, default `../external_packages/gammaLoss/tables`) |
| `config/` | example user XMLs: `jetscape_user_gamma_brick.xml` (PhotonGun + brick) and `jetscape_user_gamma_MUSIC.xml` (Trento + MUSIC) |
| `config/photon/` | run configurations for AuAu 200 GeV and PbPb 2.76 TeV (0-20%, 20-40%), pp references, brick and rate tests, and the MUSIC input file |
| `analysis/` | C++ analysis of the final state output: `AuAu200GammaAnalysis`, `PbPb2760GammaAnalysis`, `PP200GammaAnalysis`, `PP2760GammaAnalysis` |
| `data/` | experimental photon and hadron spectra and v2 (PHENIX, STAR, ALICE) used by the analysis |
| `python/` | plotting, comparison and slurm submission scripts |

Run the examples from the C-SCAPE `build` directory, for example

```bash
./runJetscape ../external_packages/gammaLoss/config/jetscape_user_gamma_brick.xml
```

The analysis programs also run from `build` and read the data from
`../external_packages/gammaLoss/data`.

## Photon status codes

All photons are written directly to the final state hadron output by
`gammaLoss`; jet photons are not passed to hadronization.

| pstat | meaning |
|-------|---------|
| 22 | hard / prompt photon from the hard process |
| 23 | photon radiated in the jet shower (e.g. MATTER) |
| 24 | thermal photon emitted by the QGP |
| 25 | thermal photon emitted by the hadron resonance gas |
| negative | absorbed in the medium, only written with `<writeAbsorbed>1` |

Hadrons from SMASH carry status 27.

## Notes on the configurations

The configurations were converted from the X-SCAPE-COMP `cparker/gamma`
branch:

- `gammaLoss` moved from `<Eloss>` into its own `<Photons>` block after the
  medium; `gammaLoss_on` and `hydro2D` were removed.
- `only_photons` was dropped since photons no longer go through
  hadronization; for photon-only runs simply leave out `<JetHadronization>`
  (configurations that used `only_photons=1` have no such block).
- Configurations without `gammaLoss` (pp references and hadron runs) have a
  `<Photons>` block with emission and reabsorption off, so hard and shower
  photons are still labeled 22 and 23.
- `min_bias` (one hard process per binary collision) is not available in
  C-SCAPE and was removed from `hydro/PbPb2760/20-40/hadrons.xml`.
- Paths into the code are relative to the C-SCAPE `build` directory. Output
  paths (`outputFilename`) and the pp baselines read by the analysis still
  point to `/scratch/user/cameron.parker/projects/runs`, and the python
  scripts assume the checkout is at
  `/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma`.
- The iSS configurations use the default `../external_packages/iSS/iSS_parameters.dat`.

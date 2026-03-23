# README for the multiscale finite element MATLAB codes

This package documents the uploaded MATLAB codes associated with the following three papers:

1. **Exponential Convergence for Multiscale Linear Elliptic PDEs via Adaptive Edge Basis Functions**  
   Yifan Chen, Thomas Y. Hou, Yixuan Wang  
   arXiv:2007.07418

2. **Exponentially convergent multiscale methods for high frequency heterogeneous Helmholtz equations**  
   Yifan Chen, Thomas Y. Hou, Yixuan Wang  
   arXiv:2105.04080

3. **Exponentially Convergent Multiscale Finite Element Method**  
   Yifan Chen, Thomas Y. Hou, Yixuan Wang  
   arXiv:2212.00823

## How the uploaded code maps to the three papers

The uploaded material is organized into **two** MATLAB archives:

- `elliptic.zip` -> code for the **elliptic** multiscale framework of arXiv:2007.07418
- `helmholtz.zip` -> code for the **Helmholtz** multiscale framework of arXiv:2105.04080

The 2022 paper arXiv:2212.00823 is a **concise review / synthesis** of the exponentially convergent multiscale finite element method (ExpMsFEM) covering both the elliptic and Helmholtz settings. There is **no separate third archive** in the upload; instead, the two codebases together serve as the computational companion to that review.

## Language and environment

The codes are written in **MATLAB** and organized as small self-contained experiment folders.

Typical requirements:

- MATLAB with sparse linear algebra support
- No PDE toolbox is required by the core scripts
- Optional: the random-coefficient generator `rando.m` uses `gmdistribution`, so regenerating the random medium may require the **Statistics and Machine Learning Toolbox**

## Recommended workflow

For most experiments:

1. Open MATLAB.
2. Change into one experiment folder.
3. Run `main.m`.

Examples:

```matlab
cd('code/periodic/Exp')
main
```

```matlab
cd('Code/case1/Exp')
main
```

Many experiments:

- compute a fine-grid reference solution with `FEM.m`
- compute a multiscale approximation with `MsFEM.m`
- report relative errors in `L2`-type and energy-type norms
- save results to `.mat` files

Some `main.m` files already save outputs, while others have the `save(...)` line commented out; if you want persistent results, simply uncomment or edit the relevant line.

## High-level organization

### 1. Elliptic code (`elliptic.zip`)

Main root inside the archive:

```text
code/
```

Main benchmark families:

- `periodic/` - periodic / oscillatory coefficient field
- `Random/` - random heterogeneous coefficient field
- `Highcontrast/` - high-contrast inclusions

Common method subfolders:

- `O(H)/` - baseline low-order coarse approximation
- `H+bubble/` - harmonic-edge plus bubble enrichment
- `Exp/` - adaptive edge-basis / exponentially convergent variant

Additional notes:

- folders with names like `copy`, `Exp copy`, `H+bubble copy` appear to be working snapshots or alternate runs
- `Nodal/` under `periodic/` is an auxiliary reference implementation rather than the main experiment entry point

### 2. Helmholtz code (`helmholtz.zip`)

Main root inside the archive:

```text
Code/
```

Main benchmark families:

- `case1/`, `case2/`, `case3/`, `case4/` - four main Helmholtz test families
- `HIgh wavenumber/` - fine-grid FEM driver for a high-wavenumber Robin-boundary test
- `Mie/` - fine-grid FEM driver for a Mie-type heterogeneous benchmark

Common method subfolders inside each case:

- `H+bubble/` - bubble-enriched multiscale basis
- `Exp/` - exponentially convergent edge-basis variant
- `H+bubble_conj/`, `Exp_conj/` - conjugate/adjoint-enriched variants used in the Helmholtz experiments
- `Nodal_can_be_regarded_as_H+bubble_when_N_e=0/` - nodal reference implementation

Additional notes:

- folders named `case1 copy`, `case2 copy`, etc. look like duplicated working snapshots
- the capitalization `HIgh wavenumber` is preserved from the original files

## Common file roles

Across both archives, the same naming pattern appears repeatedly.

- `main.m`  
  Experiment driver. Usually sets mesh parameters, wavenumber if needed, calls `FEM` and `MsFEM`, and computes relative errors.

- `FEM.m`  
  Fine-grid reference finite element solver.

- `MsFEM.m`  
  Multiscale solver for the corresponding basis family.

- `elementstiff*.m`  
  Local stiffness / load assembly on elements or coarse patches.

- `basefun.m`, `basefun1.m`  
  Local basis construction utilities.

- `harmext.m`  
  Harmonic extension component.

- `bubble.m`  
  Bubble-part construction.

- `restrict.m`  
  Restriction / edge data handling.

- `loc2glo.m` and `cellVertices.m`  
  Mesh indexing and geometry utilities.

- `afun.m`, `ffun.m`  
  Coefficient and forcing definitions.

For Helmholtz, there are additional coefficient files:

- `betafun.m` - boundary / Robin coefficient
- `gfun.m` - boundary source term
- `vfun.m` - reaction / refractive-index type coefficient
- `bc.m` - boundary indexing helper

## Reproducing the main experiments

### Elliptic experiments

Suggested entry points:

- `code/periodic/O(H)/main.m`
- `code/periodic/H+bubble/main.m`
- `code/periodic/Exp/main.m`
- `code/Random/H+bubble/main.m`
- `code/Random/Exp/main.m`
- `code/Highcontrast/H+bubble/main.m`
- `code/Highcontrast/Exp/main.m`

Typical parameters inside `main.m`:

- `N_c` - number of coarse elements per spatial direction
- `N_f` - number of fine elements per coarse element
- loop index `i = 1:7` - number of edge basis functions / local enrichment level

Typical saved outputs:

- `L` - relative `L2`-type error
- `H` - relative energy-type error

### Helmholtz experiments

Suggested entry points:

- `Code/case1/Exp/main.m`
- `Code/case1/H+bubble/main.m`
- `Code/case1/Exp_conj/main.m`
- similarly for `case2`, `case3`, `case4`

Typical parameters inside `main.m`:

- `N_c` - coarse mesh resolution
- `N_f` - fine elements per coarse element
- `k0` - wavenumber
- loop index `i = 1:7` - local edge enrichment level

Remarks:

- Several Helmholtz `main.m` files have a commented-out `save(...)` command. Uncomment it to store error arrays.
- `case4` loads coefficient data from `u.mat`, `v.mat`, and `beta.mat`.
- `HIgh wavenumber/` and `Mie/` do not have a `main.m`; instead, they expose a direct FEM solve through `FEM(N_f, k0)` after editing the coefficient files if desired.

## How to modify the PDE setup

### Elliptic archive

Edit the following files inside a chosen experiment folder:

- `afun.m` - diffusion coefficient
- `ffun.m` - right-hand side (when present)

Examples included in the archive:

- oscillatory periodic coefficients
- random coefficient fields interpolated from `u.mat`
- high-contrast inclusion patterns

### Helmholtz archive

Edit:

- `afun.m` - diffusion / material coefficient
- `vfun.m` - reaction / refractive coefficient
- `betafun.m` - Robin boundary coefficient
- `gfun.m` - boundary source
- `ffun.m` - interior source term

Some cases implement:

- homogeneous media with oscillatory boundary excitation
- heterogeneous refractive index fields
- high-contrast inclusions
- random media loaded from `.mat` arrays

## Practical notes

- Many folders are duplicated with names containing `copy`. For a clean starting point, use the folders **without** `copy` in the name.
- Several scripts create figures with `surf(...)`; if you are running in batch mode, you may want to disable plotting.
- The codes are written as experiment scripts rather than as a packaged library, so each folder is intended to be run independently.
- Relative path assumptions matter: run a script **from inside its own folder** so that local `.mat` files and helper functions are found correctly.

## Suggested citation

If you use these codes, please cite the corresponding paper(s):

- arXiv:2007.07418 for the elliptic adaptive edge-basis method
- arXiv:2105.04080 for the Helmholtz multiscale method
- arXiv:2212.00823 for the concise ExpMsFEM review

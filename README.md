# README for the elliptic multiscale finite element code

This archive contains the MATLAB code corresponding primarily to:

**Exponential Convergence for Multiscale Linear Elliptic PDEs via Adaptive Edge Basis Functions**  
Yifan Chen, Thomas Y. Hou, Yixuan Wang  
arXiv:2007.07418

It is also part of the broader ExpMsFEM framework summarized in:

**Exponentially Convergent Multiscale Finite Element Method**  
arXiv:2212.00823

## What this code does

The code implements multiscale finite element methods for second-order linear elliptic PDEs with rough coefficients. The core workflow is:

1. solve a fine-grid reference problem using `FEM.m`
2. solve a reduced multiscale problem using `MsFEM.m`
3. compare the two solutions in relative `L2`-type and energy-type norms

The experiment folders are organized by **coefficient family** and **basis construction strategy**.

## Archive layout

```text
code/
|- periodic/
|  |- O(H)/
|  |- H+bubble/
|  |- Exp/
|  \- Nodal/
|- Random/
|  |- H+bubble/
|  \- Exp/
|- Highcontrast/
|  |- H+bubble/
|  \- Exp/
\- ... additional folders named "copy", which are duplicate working snapshots
```

### Main coefficient families

- `periodic/`  
  Oscillatory periodic / multiscale coefficient benchmark.

- `Random/`  
  Random heterogeneous medium. Some folders load a precomputed coefficient field from `u.mat`. A helper `rando.m` is also included for generating random data.

- `Highcontrast/`  
  High-contrast inclusions benchmark.

### Main basis families

- `O(H)/`  
  Baseline coarse approximation.

- `H+bubble/`  
  Harmonic-edge basis enriched with the bubble part.

- `Exp/`  
  Adaptive edge-basis construction giving exponential convergence behavior.

- `Nodal/`  
  Auxiliary nodal reference code; not the main entry point used by the experiments.

## Recommended starting points

To reproduce the main reported trends, start with the folders **without** `copy` in the name:

```matlab
cd('code/periodic/Exp'); main
cd('code/periodic/H+bubble'); main
cd('code/periodic/O(H)'); main

cd('code/Random/Exp'); main
cd('code/Random/H+bubble'); main

cd('code/Highcontrast/Exp'); main
cd('code/Highcontrast/H+bubble'); main
```

## Typical experiment parameters

Most `main.m` files define:

- `N_c` - number of coarse elements per spatial direction
- `N_f` - number of fine elements per coarse element
- `i = 1:7` - enrichment level / number of edge basis functions

Example pattern:

```matlab
N_c = 32;
N_f = 32;
[result,K,C] = FEM(N_c*N_f);
for i = 1:7
    erro = result - MsFEM(N_c, i, N_f);
    ...
end
```

## Common outputs

Most drivers compute:

- `L` - relative `L2`-type error
- `H` - relative energy-type error

and save them in files such as:

- `eg1_method1_m1to7.mat`
- `eg1_method2_m1to7.mat`
- `eg1_method3_m1to7.mat`
- `eg2_method2_m1to7.mat`
- `eg2_method3_m1to7.mat`
- `eg3_method2_m1to7_*.mat`
- `eg3_method3_m1to7_*.mat`

The naming convention reflects:

- `eg1`, `eg2`, `eg3` -> benchmark family
- `method1`, `method2`, `method3` -> basis family
- `m1to7` -> enrichment level sweep

## Main files in each experiment folder

- `main.m`  
  Driver script for one benchmark/method pair.

- `FEM.m`  
  Fine-grid reference finite element solve.

- `MsFEM.m`  
  Multiscale solve for the selected basis family.

- `elementstiff.m`, `elementstiff1.m`, `elementstiff2.m`  
  Local stiffness and load assembly.

- `basefun.m`, `basefun1.m`  
  Local basis construction.

- `harmext.m`  
  Harmonic extension.

- `bubble.m`  
  Bubble-part computation.

- `restrict.m`  
  Restriction / edge-data operator.

- `loc2glo.m`, `cellVertices.m`  
  Mesh indexing and geometry helpers.

- `afun.m`  
  Diffusion coefficient.

- `ffun.m`  
  Right-hand side, when present.

## Modifying the model problem

To change the PDE setup in a folder, edit:

- `afun.m` to change the coefficient field
- `ffun.m` to change the forcing term

Examples already present in the archive include:

- highly oscillatory periodic coefficients
- random media interpolated from a sampled field
- high-contrast inclusions

If you want to regenerate the random coefficient rather than using the provided `u.mat`, inspect `rando.m`. Note that it uses `gmdistribution`, which may require the Statistics and Machine Learning Toolbox.

## Notes

- Run scripts from inside their own folders so that local helper files and `.mat` data are found correctly.
- Folders with `copy` in the name appear to be alternate snapshots or repeated runs; use them only if you specifically want those variants.
- Some folders contain plotting utilities such as `plot_a.m` and `myprint.m`.

## Citation

Please cite:

- arXiv:2007.07418 for the elliptic adaptive edge-basis method
- arXiv:2212.00823 for the broader ExpMsFEM overview

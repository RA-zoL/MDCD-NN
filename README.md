# MDCD-NN

This repository will contain validation and application data to the paper Expanding Frontiers of Complex Reaction Network Exploration through a General Reactive Machine Learning Potential, by S. Li et al.

## Folder: mdcdnn

This folder contains MDCD-NN model ensembles (4 models for `l=2` and 4 for `l=3`).

## Folder: Gaussian_interface

This folder contains the interface scripts for MDCD-NN to Gaussian 16.

To use these scripts, first, launch `SevennServer.py`:
```
Usage: SevennServer.py N model_1 model_2 ... model_N P
```
The `SevennServer.py` will be initialized with `N` models provided by `model_1` to `model_N`. After the initialization succeeded, a suggestive file `P.initialised` will be generated.

Then, MDCD-NN-based energy, gradient and hessian (numerical) calculations can be invoked by Gaussian 16 with the `SevennClient.py` script through the `external` keyword (https://gaussian.com/external/).
```
Example: #p opt=(nomicro,calcfc,recalc=10) nosymm external='SevennClient.py'
```
Finally, after all the calculations complete, one can terminate `SevennServer.py` simply by execuating `SevennTerminator.py`.

## Folder: benchmark

`RGD1_s.wB97XD4.xyz`: `gCP(TZ)-wB97X-D4/def2-TZVP` calculated elementary reaction steps selected from the RGD1 dataset. Structures are sequenced by the order of `reactant_i`, `product_i`, `TS_i`, `reactant_i+1`, ...
`Cyclo_s.wB97XD4.xyz`: `gCP(TZ)-wB97X-D4/def2-TZVP` calculated elementary reaction steps selected from the Cyclo[3+2] dataset. Structures are sequenced by the order of `reactant_A_i`, `reactant_B_i`, `product_i`, `TS_i`, `reactant_A_i+1`, ...

## Folder: validation

All files with the name `rxn_{Tag:03d}.gjf` are starting reactant structures for 100 textbook elementary reactions investigated in this work. The files are in the Gaussian input `gjf` format except several additional lines behind the cartesian coordinates to specify the reaction coordinate (RC) used in CD calculations.

For example, the file `rxn_001.gjf` contains three additional lines:
```
Refene: Auto
FixBond1: 2 15 3.19 3.04 2.89 2.74 2.58 2.43 2.28 2.13 1.98 1.82 1.67 1.52
FixBond2: 3 18 3.19 3.04 2.89 2.74 2.58 2.43 2.28 2.13 1.98 1.82 1.67 1.52
```
The first line specifies the reference energy in `a.u.` for this CD calculation, where `Auto` means to perform free optimization for the starting structure to obtain its reference energy.

The following lines starting with `FixBond{i}:` have `N+2` parameters in each line, with `N` as the number of a series of constrained optimization jobs in this CD calculation. The first two parameters indicate the atom pair `(i,j)` (starting from 1) and the resting `N` float numbers $r_1, r_2, ..., r_N$ represent that the atomic distance $d(i,j)$ will be fixed at $r_1, r_2, ..., r_N$ in `Angstrom` in each constrained optimization. If the value $r_k$ is 0.00, the distance of atom pair `(i,j)` is kept to its value in the previous constrained optimization job ($r_{k-1}$) for $k>0$ or in the starting structure for $k=0$.

In this case, we perform 12 constrained optimization jobs. The atomic distances $d(2,15)$ and $d(3,18)$ are fixed at 3.19, 3.04, 2.89, 2.74, 2.58, 2.43, 2.28, 2.13, 1.98, 1.82, 1.67, 1.52 `Angstrom`, respectively.

The tarball `MDCD-NN.tar.gz` contains the Gaussian 16 outputs of MDCD-NN computed frequency analysis `freq` jobs for each reactant, product and transition state (TS) of the 100 elementary reactions.

The tarball `DFT.tar.gz` contains the Gaussian 16 outputs of `gCP(TZ)-wB97X-D4/def2-TZVP` computed `opt freq` jobs for each reactant, product and transition state (TS) of the 100 elementary reactions.

## Folders: vindoline, endiandric and proline

All files with the suffix `cfms.xyz` are MDCD-NN-optimized conformers of intermediates and TSs. Each file represents one structural isomer, with multiple frames of structures within representing all conformational isomers of this isomer, sorted in order of relative Gibbs free energy from lowest to highest. The xyz-formatted comment line for each conformational isomer contains four numbers, from left to right, the electronic and Gibbs free energies in `a.u.`, and the relative electronic and Gibbs free energies in `kcal/mol` over the lowest-energy conformer of this structural isomer calculated with MDCD-NN, respectively.

All files with the name `DFT.xyz` are `gCP(TZ)-wB97X-D4/def2-TZVP` optimized structures for the lowest-energy conformers involved in the MDCD-NN-calculated reaction networks. The name of each structure is noted in the xyz-formatted comment line.

## Folder: benzo

`benzo_dft.xyz`: `gCP(TZ)-wB97X-D4/def2-TZVP` optimized structures.
`benzo_nnp.xyz`: MDCD-NN optimized structures involved in the double hydrogen transfer in 2,5-diamino-1,4-benzoquinone.

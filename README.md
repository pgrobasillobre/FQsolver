# FQSolver: an Open Source Code for Electrostatic Potentials, Fields, and Fluctuating Charges

<p align="center">
  <!-- FQSolver logo placeholder. Replace this path/URL when the logo is ready. -->
  <img src="docs/_static/FQSolver.png" width="600">
</p>

## Table of Contents

- [About](#about)
- [Theoretical Framework](#theoretical-framework)
- [Installation](#installation)
- [Usage](#usage)
- [Input Files](#input-files)
- [Tests](#tests)
- [License](#license)
- [Contact](#contact)

## About

**FQSolver** is a C++ code for computing electrostatic quantities from a quantum-mechanical density stored in **CUBE format**.

It currently supports three main tasks:

1. **Integration of electron density from a CUBE file**
2. **Computation of electrostatic potential and electric field at solvent sites**
3. **Computation of solvent charges through the Fluctuating Charges approach (FQ)**

Solvent geometries can be provided as:

- **XYZ files** for potential and electric-field computations
- **PDB files**, with residue/group and atom-name filters, for FQ charges, potential, and electric-field computations

FQSolver is designed for high-performance calculations using **OpenMP** and **BLAS/LAPACK** linear algebra routines.

## Theoretical Framework

All quantities in FQSolver are expressed internally in atomic units.

Given a solute electron density $\rho(\mathbf{r})$, FQSolver can compute the electrostatic potential at solvent coordinates $\mathbf{R}_i$:

$$
V(\mathbf{R}_i) =
\int d\mathbf{r} \ \frac{\rho(\mathbf{r})}{|\mathbf{R}_i - \mathbf{r}|}
$$

and the corresponding electric field:

$$
\mathbf{E}(\mathbf{R}_i) = -\nabla V(\mathbf{R}_i)
$$

For FQ calculations, solvent charges are obtained by solving a linear system based on the Fluctuating Charges model. The FQ implementation follows the approach and parametrization described in:

> Tommaso Giovannini, Alessandra Puglisi, Matteo Ambrosetti, and Chiara Cappelli,  
> *Journal of Chemical Theory and Computation* **2019**, 15 (4), 2233-2245.  
> DOI: [10.1021/acs.jctc.8b01149](https://doi.org/10.1021/acs.jctc.8b01149)

The FQ linear system has the block structure:

$$
\begin{pmatrix}
T^{qq} & C^T \\
C      & 0
\end{pmatrix}
\begin{pmatrix}
q \\
\lambda
\end{pmatrix}
=
\begin{pmatrix}
-V - \chi \\
Q
\end{pmatrix}
$$

where:

- $T^{qq}$ is the charge-charge interaction matrix
- $C$ imposes charge conservation constraints per solvent molecule
- $q$ are the fluctuating charges
- $\lambda$ are Lagrange multipliers
- $V$ is the electrostatic potential at solvent sites
- $\chi$ is the atom-type electronegativity parameter
- $Q$ is the molecular charge constraint

At present, the implemented FQ parametrization is:

- `giovannini`

and the implemented FQ interaction kernel is:

- `gaus`

## Installation

FQSolver requires:

- CMake 3.10 or higher
- C++20-compatible compiler
- BLAS/LAPACK libraries
- Python 3
- OpenMP support, recommended

To configure the project, run:

```bash
./setup.sh
```

By default, `setup.sh` tries to configure an OpenMP build first. If OpenMP is not available, it falls back to a serial build.

To require OpenMP explicitly:

```bash
./setup.sh -omp
```

To force a serial build:

```bash
./setup.sh --serial
```

To compile:

```bash
cd build
make -j
```

This generates the executable:

```bash
build/FQSolver
```

## Usage

From a test or working directory containing an input file:

```bash
/path/to/FQSolver input_file.inp
```

To select the number of OpenMP threads:

```bash
/path/to/FQSolver input_file.inp -omp 8
```

Output files are written to log and text files. Text results are stored in:

```text
FQSolver_results/
```

## Input Files

### Potential from Density

```text
what: potential
solvent: solvent/water_solvent.xyz
density: density/zn-pc_scf_dens.cube
cutoff: 1e-03
```

### Field from Density

```text
what: field
solvent: solvent/water_solvent.xyz
density: density/zn-pc_scf_dens.cube
cutoff: 1e-03
```

### Potential and Field from Density

```text
what: potential+field
solvent: solvent/water_solvent.xyz
density: density/zn-pc_scf_dens.cube
cutoff: 1e-03
```

### PDB Solvent Parsing

For PDB files, solvent atoms can be filtered by group and atom names:

```text
what: potential+field
solvent: solvent/frame186_18.pdb
group: sol x
read atoms: OW, HW1, HW2
density: density/pdi_trans1.cube
cutoff: 1e-03
```

Here:

- `group` selects the residue/group string to read
- `read atoms` selects the atom labels to keep

### FQ Charges from Density

```text
what: fq
parametrization: giovannini
kernel: gaus

solvent: solvent/frame186_18.pdb
group: sol x
read atoms: OW, HW1, HW2

density: density/pdi_trans1.cube
cutoff: 1e-03

debug: 0
```

The FQ charge result file contains:

```text
Atom, Mol. Index, X (Å), Y (Å), Z (Å), Potential (a.u.), Charge (a.u.)
```

### Density Integration

```text
integrate cube file: zn-pc_scf_dens.cube
```

## Tests

After running `./setup.sh`, a parallel test runner is generated in the build directory.

To run the full test suite:

```bash
cd build
./run_tests.sh
```

The script runs CTest using all detected processors by default.

You can also run CTest manually:

```bash
ctest --output-on-failure
```

or in parallel:

```bash
ctest -j 8 --output-on-failure
```

Example input files and references are located in:

```text
tests/
```

## License

FQSolver is licensed under the **GNU General Public License v3.0**.

## Contact

For issues or contributions:

- Email: **pgrobasillobre@gmail.com**
- GitHub issues: https://github.com/pgrobasillobre/FQsolver/issues

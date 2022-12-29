# dmrg_CPP

Archived programs for calculating one-dimensional models about condensed matter physics using density matrix renormalization group (DMRG) and exact diagonalization.
This repository will be archived.

## How to Use
First, you go enter `dmrg_CPP/lib` and do make

```bash
cd lib
make
```

If you use Linux systems, please comment out Mac OS clang commands
https://github.com/K-Suzuki-Jij/dmrg_CPP/blob/9830c86fbf079d3416be74aa4f6893ec3c4d181a/lib/makefile#L12-L14

and then uncomment Linux command for using g++
https://github.com/K-Suzuki-Jij/dmrg_CPP/blob/9830c86fbf079d3416be74aa4f6893ec3c4d181a/lib/makefile#L16-L18

and for Intel oneAPI
https://github.com/K-Suzuki-Jij/dmrg_CPP/blob/9830c86fbf079d3416be74aa4f6893ec3c4d181a/lib/makefile#L20-L22

We check the programes work on Ubuntu 20.04.5 LTS (GNU/Linux 5.15.0-56-generic x86_64)

Next, you go enter any directries under `dmrg_CPP/main/dmrg`. Here, we go enter `dmrg_CPP/main/dmrg/EKLM` and do make

```bash
cd dmrg_CPP/main/dmrg/EKLM
make
```

Again here, If you use Linux systems, please comment out Mac OS clang commands.
https://github.com/K-Suzuki-Jij/dmrg_CPP/blob/9830c86fbf079d3416be74aa4f6893ec3c4d181a/main/dmrg/EKLM/makefile#L10-L13
and then uncomment Linux command.
Finally, you can execute program appeared in `dmrg_CPP/main/dmrg/EKLM/build/1D_EKLM_DMRG.out`

```
./build/1D_EKLM_DMRG.out
```

Calculation log will be appeared.
<img width="1365" alt="Screenshot 2022-12-30 at 0 30 49" src="https://user-images.githubusercontent.com/78338408/209974953-a4984df4-723d-48de-86dd-287ca8145e3d.png">

Results will be written in `dmrg_CPP/main/dmrg/EKLM/result` and status will be written in `dmrg_CPP/main/dmrg/EKLM/SML_out`

## Programs by DMRG
* dmrg_CPP/main/dmrg/AKLM: Anisotropic Kondo lattice model. (under construction)
* dmrg_CPP/main/dmrg/AKLM_TVF: Anisotropic Kondo lattice model under the transeverse fields. (under construction)
* dmrg_CPP/main/dmrg/EKLM: Extended Kondo lattice model.
* dmrg_CPP/main/dmrg/HUBBARD: Hubbard model. (under construction)
* dmrg_CPP/main/dmrg/TAKLM: Two-channel anisotropic Kondo lattice model. (under construction)
* dmrg_CPP/main/dmrg/XXZ: XXZ model. (under construction)

## Programs by Exact Diagonalization
* dmrg_C/main/exact_diag/XXZ: XXZ model.

## Environment
* Mac OS 13.1 on Apple M1 Max
* Mac OS 13.1 on Intel CPU
* Ubuntu 20.04.5 LTS on Ryzen 7950x

# Overview

`galotfa`: galactic on-the-fly analysis, a library for on-the-fly analysis of
galactic simulations.

At present, mainly focused on disk or bar galaxies. This library is implemented
in C++ and offers C API. It uses `HDF5` for data storage, and `MPI` for parallelization.

---

## Dependencies

- A C++ compiler that supports C++17.
- [`HDF5`](https://www.hdfgroup.org/solutions/hdf5/)
- [`gsl`](https://www.gnu.org/software/gsl/)
- A `MPI` library (e.g. [`OpenMPI`](https://www.open-mpi.org/)).
- `cmake` >= 3.12

- [`toml++`](https://marzer.github.io/tomlplusplus/#mainpage-example) for `toml`
  parser, open source with MIT license, already included in this project.

---

## Installation

---

## Usage

---

## Examples

The structure of the HDF5 log file:

- Orbit:
  Restore several orbits information:
  - Times: m by 1 (m: number of log steps)
  - Masses: m by n by 1 (n: number of orbits)
  - ParticleIDs: m by n by 1
  - Coordinates: m by n by 3
  - Velocities: m by n by 3

---

## Future Work

---

## Acknowledgements

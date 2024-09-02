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
  Restore several orbits as datasets, each names as "Particle-X", where X is the id of the particle.
  Each dataset has N x 7 numbers: N is the number of logged steps, 7 for 1 time, 3 coordinates, and
  3 velocities in Cartesian coordinate frame.

---

## Future Work

---

## Acknowledgements

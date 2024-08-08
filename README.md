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
- [`xmake`](https://xmake.io/#/about/introduction) for building, it's an awesome
  modern and fast build system. The installation is very easy:

```bash
  curl -fsSL https://xmake.io/shget.text | bash
```

Note that if you want to run the unit tests with `xmake`, you need `xmake`
later than 2.8.4 !

- [`toml++`](https://marzer.github.io/tomlplusplus/#mainpage-example) for `toml`
  parser, open source with MIT license, already included in this project.

---

## Installation

---

## Usage

---

## Examples

---

## Future Work

---

## Acknowledgements

# Nevanlinna

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/SpM-lab/Nevanlinna.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/SpM-lab/Nevanlinna.jl/dev)
[![Build Status](https://github.com/SpM-lab/Nevanlinna.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SpM-lab/Nevanlinna.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Installation
The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Nevanlinna
```

This will install a Command Line Interface (CLI) script, `nevanlinna`, into `~/.julia/bin`.
You can add this directory to your PATH in a bash shell by adding the following line into `~/.bashrc`:

```bash
export PATH="$HOME/.julia/bin:$PATH"
```

This command needs input parameter TOML file.
These files can be downloaded from [here](https://github.com/SpM-lab/Nevanlinna.jl/blob/comonicon/comonicon/bare/config.toml).
TO DO:the link must be modified after marged to main branch!


### How to run examples

The notebooks in which reproduce examples in our paper are located in notebooks derectory.
These examples includes:  
- $\delta$-function
- Gaussian
- Lorentzian
- Two peak
- Kondo resonance
- tractable Hubbard gap
- challenging Hubbard gap

- compare 64-bit and 128-bit
- Hamburger moment problem

To run our code, please ensure that the following packages is already installed to Julia
- Nevanlinna
- Plots
- LaTexStrings
- SparseIR

### Manual installation from source (advanced)

You should almost never have to do this, but it is possible to install Nevanlinna.jl from source as follows:
```sh
git clone https://github.com/nogaki/Nevanlinna.jl.git
julia -e "import Pkg; Pkg.add(path=\"Nevanlinna.jl\")"
```
This is *not* recommended, as you will get the unstable development version and no updates.

###To install CLI
```
julia --project deps/build.jl install
```

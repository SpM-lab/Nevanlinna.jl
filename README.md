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


### How to run examples

TO BE WRITTEN

### Manual installation from source (advanced)

You should almost never have to do this, but it is possible to install Nevanlinna.jl from source as follows:
```sh
git clone https://github.com/nogaki/Nevanlinna.jl.git
julia -e "import Pkg; Pkg.add(path=\"Nevanlinna.jl\")"
```
This is *not* recommended, as you will get the unstable development version and no updates.

#To install CLI
```
julia --project deps/build.jl install
```

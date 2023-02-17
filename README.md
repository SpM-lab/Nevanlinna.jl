# Nevanlinna

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/SpM-lab/Nevanlinna.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://github.com/SpM-lab/Nevanlinna.jl/dev)
[![Build Status](https://github.com/SpM-lab/Nevanlinna.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SpM-lab/Nevanlinna.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Installation
The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Nevanlinna
```

### Manual installation from source

You should almost never have to do this, but it is possible to install Nevanlinna.jl from source as follows:
```sh
git clone https://github.com/nogaki/Nevanlinna.jl.git
julia -e "import Pkg; Pkg.add(path=\"Nevanlinna.jl\")"
```
This is *not* recommended, as you will get the unstable development version and no updates.

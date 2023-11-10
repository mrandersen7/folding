<div align="center">

# **Protien Folding Subnet** <!-- omit in toc -->
[![Discord Chat](https://img.shields.io/discord/308323056592486420.svg)](https://discord.gg/bittensor)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 

---

## The Incentivized Internet <!-- omit in toc -->

[Discord](https://discord.gg/bittensor) • [Network](https://taostats.io/) • [Research](https://bittensor.com/whitepaper)
</div>

---
- [Introduction](#introduction)
- [Installation](#installation)
  - [Before you proceed](#before-you-proceed)
- [Background: What is protein folding?](#background)
- [Description](#description)
- [Features](#features)

- [License](#license)



## Introduction

  This is the start of the protein folding subnet. This subnet uses the Gromacs software to simulate molecular dynamics of proteins. We take a known initial 3D structure, and put in a cell-like environment and simluate it to know its end form. This is an essential step in the protein folding process and an entry point to many other high level techniques.

General informaiton about gromacs can be found here: https://manual.gromacs.org/2023.2/index.html


## Installation
GROMACS
You will need two packages to run this miner. Gromacs itself, and then a gromacs wrapper to make the base functions more python friendly. You can find the install process and requirements for the latest version of gromacs here:
- https://manual.gromacs.org/2023.2/install-guide/index.html

However, I found package managers make the process much simpler based on your preffered workflow:
- Conda install: conda install -c conda-forge gromacs
- Brew install: brew install gromacs

GromacsWrapper
For the most part, this leaves the base syntax intact. Installation instructions and more can be found here: https://gromacswrapper.readthedocs.io/en/latest/installation.html
- Conda install: conda install -c conda-forge gromacswrapper
- pip install: pip install GromacsWrapper


### Before you proceed

  Comlpexity of the problem aside, one of the current barriers to protein folding is computational ability. The proceess involved are complex and take time even with state of the art systems. In contrast to other subnets, this specific miner process can take hours. There are 3 mdrun functions, each will output a projected time when run. The first two will be about the same length, with the third taking about an order of magnitude longer. However, this is to run each run to 50,000 steps. If you wish to run a shorter simluation, you can use the maxh argument in each mdrun to limit how long you would like to spend (You can also limit steps if you prefer). Your simulations can be resumed later, or even continued with the "incomplete" data, but you should be aware of the timescale involved with undertaking this process first
  
  
## Background  
  
  Proteins are the biological molecules that "do" things, they are the molecular machines of biochemistry. Enzymes that break down food, hemoglobin that carries oxygen in blood, and actin filaments that make muscles contract are all proteins. They are made from long chains of amino acids, and the sequence of these chains is the information that is stored in DNA. However, its a large step to go from a 2D chain of amino acids to a 3D structure capable of working. 

  The process of this 2D structure folding on itself into a stable, 3D shape in a cell is called protein folding. For the most part, this process happens naturally and the end structure is in a much lower free energy state than the string. Like a bag of legos though, its not enough to just know the building blocks being used, its the way they're supposed to be put together that matters. "Form defines function" is a common phrase in biochemsitry, and it is the quest to determine form, and thus function of proetins, that makes this process so important to understand and simulate. 

  Understanding how specific proteins fold unlocks the ability to cure many ailments. Folding@Home, a distribtred computing community dedicated to simulating protein folding, was able to help design a treatment for SARS-covid-19 by identifying a unique folding patter in the spike protein of the virus that left it open to interference. Understading how beta amyloid plaques fold, and this misfold, is essential to understanding how Alzheimers Disease develops and potential treamtent points.


## Description

Validators receive an input of a pdb file for a protein, with optional inputs for force field, simluation box shape, and molecular dynamics files (Defaults are Charmm27, a rhombic dodecahedron, and provided charmm27 files). This is pre-processed, with the results sent to the miner

The miner then runs 3 simluations, first stabilizing for temperature, then pressure, then running the final "production" run. These files are then sent to the validator for evaluation

The validator is going to evaluate the files for protein stablity to ensure things were run correctly, then rank outputs based on the lowest free energy state. Rewards will be distributed accordingly to all miners who particpated in the query and delivered valid results



## Features

- Gromacs itself is a rather robust package. There are specific guides and functions if you wish to parallelize your processing or run these computations off of a GPU to speed things up





## License

[INSERT GROMACS LICENSING INFORMATION HERE]



This repository is licensed under the MIT License.
```text
# The MIT License (MIT)
# Copyright © 2023 Yuma Rao

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the “Software”), to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of
# the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
```

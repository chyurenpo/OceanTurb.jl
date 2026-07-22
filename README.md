# OceanTurb.jl

| **Documentation**             | **Build Status**                    | **License** |
|:-----------------------------:|:-----------------------------------:|:-----------:|
| [![docs][docs-img]][docs-url] | [![travis][travis-img]][travis-url] |[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://mit-license.org/)|


# OceanTurb.jl: Additional Turbulence Closure Schemes

This repository is an experimental research fork of
[OceanTurb.jl](https://github.com/glwagner/OceanTurb.jl), a Julia package
for studying turbulence models and parameterizations in ocean boundary layers.

The purpose of this fork is to implement and evaluate additional turbulence
closure schemes that are not included in the original OceanTurb.jl repository.
Current development focuses on the Mellor–Yamada Level 2.5 (MY2.5), k–ω,
and k–ε closure schemes for one-dimensional ocean boundary-layer simulations.
This fork also focuses on the development and evaluation of bottom
boundary-layer parameterizations formulated in slope-following coordinate
systems.

## Current additions and development status

The following table summarizes the current implementation status and planned
development of this research fork.

| Component | Description | Status |
|---|---|---|
| Mellor–Yamada Level 2.5 (MY2.5) closure | Implementation of a prognostic turbulent kinetic energy closure for one-dimensional ocean boundary-layer simulations | Experimental implementation completed|
| k–ω closure | Two-equation turbulence closure based on turbulent kinetic energy and specific dissipation rate | Planned |
| k–ε closure | Two-equation turbulence closure based on turbulent kinetic energy and dissipation rate | Planned |
| Sloping-coordinate formulation | Extension of the one-dimensional governing equations to slope-following coordinates | In progress |
| Sloping bottom boundary-layer parameterizations | Application and evaluation of turbulence closures for ocean bottom boundary layers over sloping topography | In progress |

## Test Case
Coming soon!


## AI-assisted development

Generative AI tools, including ChatGPT by OpenAI, were used to assist with
parts of code drafting and debugging.

The governing equations, numerical formulation, implementation, test cases,
and scientific interpretation were reviewed by the repository author, who
takes responsibility for the released code and results.



[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-url]: https://glwagner.github.io/OceanTurb.jl/dev/

[travis-img]: https://travis-ci.com/glwagner/OceanTurb.jl.svg?branch=master
[travis-url]: https://travis-ci.com/github/glwagner/OceanTurb.jl

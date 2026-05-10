# Welcome Coronograph package documentation

## About

A paraxial coronagraph model for calculating axial distances of optical 
components in Lyot configuration with external occulter.

## Installation

```bash
pip install git+https://github.com/tomas38/coronagraph.git
```

## Quick Start

```python
from coronagraph import Coronagraph

c = Coronagraph(Ra=0.05, theta_v0=0.1, theta_v1=0.2, theta_m=0.3)
```

## Contents

- Usage
  - [Coronagraph theory](coronagraph_theory.md)
  - [Vignetting](Vignetting.ipynb)

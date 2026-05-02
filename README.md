# Coronagraph

A paraxial coronagraph model for calculating axial distances of optical 
components in Lyot configuration with external occulter.

![alt text](<docs/img/Snímek obrazovky 2026-04-19 162838.png>)
![alt text](<docs/img/Snímek obrazovky 2026-04-19 162847.png>)

## Installation

```bash
pip install git+https://github.com/tomas38/coronagraph.git
```

## Quick Start

```python
from coronagraph import Coronagraph

c = Coronagraph(Ra=0.05, theta_v0=0.1, theta_v1=0.2, theta_m=0.3)
```

## Documentation

Full documentation available at: https://tomas38.github.io/coronagraph

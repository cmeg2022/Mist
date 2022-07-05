# MiSt (Microstructure Statistics) Python library

The MiSt library can be used for the analysis of both experimental and simulated microstructures. The algorithms implemented in the library have been described in the paper available at the link here.

## Installation Guide

Clone this git repository in the folder where you wish to analyze images with:

```
git clone https://github.com/cmeg2022/Mist.git
```

To install the base libraries:

```
pip install -r requirements.txt
```

To access the scripts involed : use the following import statements.

```
from Scripts import SpatialCorrelations as corr
from Scripts import plots
from Scripts import hoshenKopleman as hosh
from Scripts import convexity as conv
```

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg


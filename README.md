# gdr2_completeness
This package helps with

[1] querying GaiaDR2 TAP services in chunks of healpixel

[2] assembling approximate completeness maps per healpixel and G magnitude bin

![](data/ag_completeness.png)
For example for GaiaDR2 extinction extimates were delivered down to G=17. This is an all-sky map in Galactic coordinates showing the completeness of all sources with extinction estimate compared to all sources down to 17th Gmag.

[3] Downloading data to assess the magnitude limits per healpix

[4] Calculating the magnitude limits of Gaia per healpix for different bands

## Getting started
The jupyter [tutorial](https://github.com/jan-rybizki/gdr2_completeness/tree/master/tutorials) illustrates the basic usage of this software. It can be inspected in the github repository or you can run it interactively on your local machine.


## Installation

```
pip install git+https://github.com/jan-rybizki/gdr2_completeness.git
```
it has a few dependencies (healpy, pyvo) and was only tested using python 3


## Authors
- Jan Rybizki (MPIA, rybizki@mpia.de)

## Collaborators
- Ronald Drimmel (INAF)
- Josep Manel Carrasco (University of Barcelona)

## Attribution
A paper is in preparation, until then an [ASCL-cite](https://ascl.net/code/v/1981) is available for citation.

## License
The package is licensed under the Creative Commons [Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/).

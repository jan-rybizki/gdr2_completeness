# gdr2_completeness
This package helps with

[[1]](https://github.com/jan-rybizki/gdr2_completeness/blob/master/tutorials/%5B1%5DQuery_GDR2_in_hpx_chunks.ipynb) querying GaiaDR2 TAP services in chunks of healpixel

[[2]](https://github.com/jan-rybizki/gdr2_completeness/blob/master/tutorials/%5B2%5DCompleteness%20tutorial_gdr2_light.ipynb) assembling approximate completeness maps per healpixel and G magnitude bin

![](data/ag_completeness.png)

For example for GaiaDR2 extinction extimates were delivered down to G=17. This is an all-sky map in Galactic coordinates showing the completeness of all sources with extinction estimate compared to all sources down to 17th Gmag.

## update in version 1.0 

[[3]](https://github.com/jan-rybizki/gdr2_completeness/blob/master/tutorials/%5B3%5DDownload%20for%20Gaia%20magnitude%20limit%20over%20the%20sky.ipynb) Downloading data to assess the magnitude limits per healpix

[[4]](https://github.com/jan-rybizki/gdr2_completeness/blob/master/tutorials/%5B4%5D%20GaiaDR2%20magnitude%20limit%20per%20healpix%20in%20different%20photometric%20bands.ipynb) Calculating the magnitude limits of Gaia per healpix for different bands

![](data/g_gdr2.png)

The all-sky G magnitude limit of GaiaDR2 sources that have BP and RP color. These maps are saved in 'data/maglim_X_in_nside128.npy' for G, BP and RP magnitude limits for healpix level 7 ordered by healpix. Or you can generate them yourself with notebook [3] and [4].

## update in version 2.0

[[5]](https://github.com/jan-rybizki/gdr2_completeness/blob/master/tutorials/%5B5%5D%20RVS_selection_function.ipynb) Querying the GDR2 RVS selection function

[[6]](https://github.com/jan-rybizki/gdr2_completeness/blob/master/tutorials/%5B6%5D%20Second%20telescope%20source%20density%20contribution.ipynb) Estimating the average number of sources from the second telescope

![](data/rvs_completeness_function.png)

The all-sky completeness of the GDR2_RVS sample with respect to all sources in GDR2 for G = 12.5 mag and GRP = 12.0 mag (corresponding to GRVS = 11.8 mag) as returned by our python function.

For more details on the RVS sample selection function please have a look at our [paper](https://ui.adsabs.harvard.edu/abs/2021MNRAS.500..397R/abstract) and you can also check out [my GalaxyCoffee talk slides](https://github.com/jan-rybizki/gdr2_completeness/blob/master/data/GC_RVS_selection_function_Rybizki.pdf).

## Getting started
The jupyter [tutorials](https://github.com/jan-rybizki/gdr2_completeness/tree/master/tutorials) illustrate the basic usage of this software. It can be inspected in the github repository or you can run it interactively on your local machine.


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
If you use sofwares in the example tutorials [1] & [2] (first release) please cite [ASCL-cite](https://ascl.net/code/v/1981).
If you use software in the example tutorials [3] & [4] (v1.0) please cite [GeDR3mock paper](https://ui.adsabs.harvard.edu/abs/2020PASP..132g4501R/abstract).
If you use software in the example tutorials [5] & [6] (v2.0) please cite [GDR2 selection function](https://ui.adsabs.harvard.edu/abs/2021MNRAS.500..397R/abstract).


## License
The package is licensed under the Creative Commons [Attribution-ShareAlike 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/).

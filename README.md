# AICON

AICON is a program aims to achieve fast and accurate estimation of transport properties, such as electrical conductivity and thermal conductivity. The first version of AICON is specific for lattice thermal conductivity calculation. Based on the [Debye-Callaway model](https://link.aps.org/doi/10.1103/PhysRev.113.1046) for thermal conductivity calculation, we expand this theory and optimize its calculating process. All of these efforts makes this software be able to calculating lattice thermal conductivity quickly and accurately without needing any transcendental parameter. This is especially important for high-throughput calculation. Beside lattice thermal conductivity, one can also obtain other useful information such as contribution rate of acoustic branch and optic branch repectively, scattering rate of different scattering mechanism like normal process and umklapp process. For more information, check out our [article]().

AICON has a DOI:10.5281/zenodo.3371097, you can cite this code like this:

    Tao Fan. (2019, August 19). Baijianlu/mDCThermalC: mDCThermalC (Version v1.1.1). Zenodo. http://doi.org/10.5281/zenodo.3371097

## prerequisites
AICON is a Python module. AICON's runtime requirements are Python version 3.5 or higher, and the Python libraries [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/), [spglib](https://atztogo.github.io/spglib/) and [pymatgen](http://pymatgen.org/index.html).  All of them can be easily obtained from the [Python Package Index](https://pypi.python.org/pypi) (PyPI), using tools such as pip. They may also be bundled with Python distributions aimed at scientists, like [Anaconda](https://anaconda.org/), and with a number of Linux distributions. Here we recommend to use Anaconda so that dependencies should be resolved automatically.

## Compiling and install AICON
Users installing from source must install the dependencies first and then run:

    $ python setup.py install
    
## Running the tests
The distribution includes four examples: diamond, Si, SnSe, Mg2Si. The former three examples are also described in the AICON [paper](). Read /doc/UserManual to learn how to use this software and more information about the output files. 

# mvmdist
This package provides a class-based interface, similar to MATLAB's build-in functions for handling [Gaussian mixture models](https://de.mathworks.com/help/stats/gmdistribution.html). For further information on circular probability distributions and von Mises mixture models in particular, these papers give a comprehensive overview on the general concepts:

* [J. Bentley (2006), "Modelling Circular Data using a Mixture of von Mises and Uniform Distributions"](https://www.stat.sfu.ca/content/dam/sfu/stat/alumnitheses/MiscellaniousTheses/Bentley-2006.pdf)
* [I. S. Dhillon and S. Sra (2003), "Modeling Data using Directional Distributions"](http://www.cs.utexas.edu/users/inderjit/public_papers/tr03-06.pdf)

Installation
------------

1. Clone or download this repository and add it to your MATLAB search path.
2. Compile the *.mex-file located at `@VonMisesMixture/private` by browsing to that directory in MATLAB and running ```mex sampleVonMises.c``` from the command line.

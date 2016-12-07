# mvmdist
A Matlab package for probabilistic modeling of circular data with mixtures of von Mises distributions.

This package provides a class-based interface, similar to MATLAB's build-in functions for handling [Gaussian mixture models](https://de.mathworks.com/help/stats/gmdistribution.html). For further information on circular probability distributions and von Mises mixture models in particular, these papers give a comprehensive overview on the general concepts:

* [J. Bentley (2006), "Modelling Circular Data using a Mixture of von Mises and Uniform Distributions"](https://www.stat.sfu.ca/content/dam/sfu/stat/alumnitheses/MiscellaniousTheses/Bentley-2006.pdf)
* [I. S. Dhillon and S. Sra (2003), "Modeling Data using Directional Distributions"](http://www.cs.utexas.edu/users/inderjit/public_papers/tr03-06.pdf)

## Installation

1. Clone or download this repository and add it to your MATLAB search path.
2. Compile the *.mex-file located at `@VonMisesMixture/private` by browsing to that directory in MATLAB and running ```mex sampleVonMises.c``` from the command line.

## Examples

### Constructing a von Mises mixture model

The following code creates a von Mises mixture model with two mixture components.
```matlab
p = [0.5; 0.5];         % Mixture weights.
mu = [-pi/2; pi/2];     % Component means.
kappa = [5; 10];        % Concentration parameters of components.

vmm = VonMisesMixture(p, mu, kappa);
```

### Plotting the probability density function

Given a vector of angular values, the probability density function of a von Mises mixture model can be computed with the ```pdf()``` method provided by the ```VonMisesMixture``` class. The following code example plots the probability density function of a ```VonMisesMixture``` class instance.
```matlab
angles = linspace(-pi, pi, 1000)';  % The pdf() function expects a column-vector as input.
likelihoods = vmm.pdf(angles);
plot(angles, likelihoods); grid on;
```

### Drawing samples from a von Mises mixture model

This implementation uses the method introduced by [Barabesi (2005)](http://sa-ijas.stat.unipd.it/sites/sa-ijas.stat.unipd.it/files/417-426.pdf) to generate samples from a von Mises distribution. To speed up the sampling process, a mex-function is used by default which is significantly faster than the plain MATLAB implementation (especially for a large number of samples).
```matlab
nSamples = 10000;
samples = vmm.random(nSamples);     % Generate samples using *.mex-function.
```
If the plain MATLAB sampling method should be used, a second argument has to be passed to the ```random()``` function which sets the internal ```useMex``` flag to ```false```:
```matlab
samples = vmm.random(nSamples, false); % Generate samples using MATLAB implementation.
```

### Maximum likelihood parameter estimation

Parameter estimation is conducted via an Expectation-Maximization (EM) approach described in [Hung et al. (2012)](http://www.tandfonline.com/doi/abs/10.1080/02664763.2012.706268). The following example shows how a model can be fitted on a set of data samples drawn from an initial von Mises mixture distribution with three components.
```matlab
p = [0.3; 0.4; 0.3];         % Mixture weights.
mu = [-pi/2; 0; pi/3];       % Component means.
kappa = [5; 10; 2.5];        % Concentration parameters of components.

vmm = VonMisesMixture(p, mu, kappa);  % Initialize model.
samples = vmm.random(10000);          % Draw 10000 samples.

% Fit new model on data samples assuming 3 mixture components.
nComponents = 3;
fittedVmm = fitmvmdist(samples, nComponents);

% Plot initial and fitted distributions.
angles = linspace(-pi, pi, 1000)';
likelihoodsInitial = vmm.pdf(angles);
likelihoodsFitted = fittedVmm.pdf(angles);

plot(angles, likelihoodsInitial); hold on;
plot(angles, likelihoodsFitted); grid on;
axis([-pi, pi, 0, 1]);
```

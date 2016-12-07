function obj = fit(obj, angles, nComponents, varargin)
% FIT This function estimates the parameters of a mixture of
%   von Mises distributions using an Expectation-Maximization
%   scheme. The algorithm used in this function is a highly
%   modified version of the basic parameter estimation
%   procedure described in [1]. This algorithm uses a different
%   Maximum-Likelihood estimator for the concentration
%   parameter and computes soft assignments of the data points
%   to the mixture components.
%
% REQUIRED INPUTS:
%   angles - Nx1 vector, containing N angular values, ranged
%       between -pi and pi.
%   nComponents - Number of mixture components that should be
%       estimated.
%
% PARAMETERS:
%   ['MaxIter', maxIter] - Maximum number of iterations the
%       EM-algorithm should run (default = 100).
%   ['ErrorThreshold', errorThreshold] - Minimum error that
%       should be used as a stopping-criterion for the
%       EM-algorithm during convergence testing (default =
%       1E-4).
%
% LITERATURE:
%   [1] Hung et al. (2012): "Self-updating clustering algorithm
%       for estimating the parameters in mixtures of von Mises
%       distributions"
%
% AUTHOR:
%   Copyright (c) 2016      Christopher Schymura
%                           Cognitive Signal Processing Group
%                           Ruhr-Universitaet Bochum
%                           Universitaetsstr. 150
%                           44801 Bochum, Germany
%                           E-Mail: christopher.schymura@rub.de

% Check inputs
p = inputParser();
defaultMaxIter = 100;
defaultErrorThreshold = 1E-4;

p.addRequired('Angles', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'vector', '>=', -pi, '<=', pi}) ...
  );

p.addRequired('NComponents', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'integer', 'scalar', 'positive'}) ...
  );

p.addParameter('MaxIter', ...
  defaultMaxIter, ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'integer', 'scalar', 'positive'}) ...
  );

p.addParameter('ErrorThreshold', ...
  defaultErrorThreshold, ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'scalar', 'nonnegative'}) ...
  );

p.parse(angles, nComponents, varargin{:});

% Get number of samples
nSamples = length(p.Results.Angles);

% Check for fixed mean parameters and run conventional k-means to get an 
% initial clustering
cIdx = obj.ckmeans(p.Results.Angles, p.Results.NComponents, ...
  'Replicates', 1);

% Generate initial gamma matrix (hard assignments from k-means results)
gamma =  zeros(length(p.Results.Angles), p.Results.NComponents);
for sampleIdx = 1 : nSamples
  gamma(sampleIdx, cIdx(sampleIdx)) = 1;
end

% Normalize
gamma = gamma + eps;
gamma = bsxfun(@rdivide, gamma, sum(gamma, 2));

% Get initial parameter estimates
[muHat, kappaHat, cPropHat] = obj.estimateParameters( p.Results.Angles, gamma );

% Initialize log-likelihood and status parameters
logLik = -realmax;
converged = false;

for k = 1 : p.Results.MaxIter
  % E-Step: Update gamma-matrix
  gamma = ...
    obj.computeGamma(p.Results.Angles, cPropHat, muHat, kappaHat);
  
  % M-Step: Re-estimate the distribution parameters
  [muHat, kappaHat, cPropHat] = obj.estimateParameters(angles, gamma);
  
  % Evaluate the log-likelihood
  logLikNew = ...
    sum(log(obj.computeLikelihood(p.Results.Angles, cPropHat, muHat, ...
    kappaHat)));
  
  % Check for convergence
  if abs(logLik - logLikNew) < p.Results.ErrorThreshold
    % If converged, save parameters
    converged = true;
    
    % Terminate EM-algorithm
    break;
  else
    % If not converged, update log-likelihood and proceed
    % with next iteration.
    logLik = logLikNew;
  end
end

% Initialize model
obj = VonMisesMixture(cPropHat, muHat, kappaHat);
obj.logLikelihood = logLik;
obj.hasConverged = logical(converged);
obj.nIterations = k;

end

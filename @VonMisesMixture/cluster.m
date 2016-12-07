function [idx, nLogLik, posteriors] = cluster(obj, angles)
% CLUSTER This functions partitions a set of N angles into K
%   clusters, where K is determined by the number of mixture
%   components of the von Mises mixture distribution defined by
%   obj.
%
% REQUIRED INPUTS:
%   angles  - Nx1 vector, containing N angular values, ranged
%       between -pi and pi.
%
% OUTPUTS:
%   idx           - Nx1 vector, containing the cluster indices for each
%                   observation. The cluster index gives the component with
%                   the largest posterior probability for the observation,
%                   weighted by the component probability.
%   nLogLik       - The negative log-likelihood of the data.
%   posteriors    - NxK matrix, representing the posterior
%                   probabilities of each component for each observation.
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

p.addRequired('angles', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'column', '>=', -pi, '<=', pi}));

p.parse(angles);

% Compute the negative log-likelihood
nLogLik = -sum(log(obj.computeLikelihood(p.Results.angles, ...
  obj.componentProportion, obj.mu, obj.kappa)));

% Compute posterior probabilities
posteriors = obj.computeGamma(p.Results.angles, ...
  obj.componentProportion, obj.mu, obj.kappa);

% Compute cluster indices
[~, idx] = max(posteriors, [], 2);

end

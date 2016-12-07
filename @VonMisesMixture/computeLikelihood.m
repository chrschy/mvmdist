function [likelihood, compLik] = computeLikelihood(angles, cProp, mu, kappa)
% COMPUTELIKELIHOOD This function returns the total likelihood value and
%   the likelihood for each mixture component for a given set of angles.
%
% REQUIRED INPUTS:
%   angles          - Nx1 vector, containing N angular values, ranged
%                     between -pi and pi.
%   cProp           - Kx1 vector, representing the proportions of all
%                     mixture components in the distribution.
%   mu              - Kx1 vector, containing the circular means of all K
%                     components of the mixture distribution. The range of
%                     each element is bounded in [-pi, pi].
%   kappa           - Kx1 vector, containing the concentration parameters
%                     of the mixture distribution.  All elements are
%                     nonnegative.
%
% OUTPUTS:
%   likelihood      - Nx1 vector containing the likelihood values
%                     for all input angles.
%   compLik         - NxK matrix, containing the component-wise
%                     likelihood values for all input angles.
%
% AUTHOR:
%   Copyright (c) 2016      Christopher Schymura
%                           Cognitive Signal Processing Group
%                           Ruhr-Universitaet Bochum
%                           Universitaetsstr. 150
%                           44801 Bochum, Germany
%                           E-Mail: christopher.schymura@rub.de

p = inputParser();

p.addRequired('Angles', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'column', '>=', -pi, '<=', pi}) ...
  );

p.addRequired('CProp', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'vector'}) ...
  );

p.addRequired('Mu', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'vector','>=', -pi, '<=', pi, 'numel', length(cProp)}) ...
  );

p.addRequired('Kappa', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'vector', 'nonnegative', 'numel', length(cProp)}) ...
  );

p.parse(angles, cProp, mu, kappa);

% Get total number of samples and mixtures and initialize component-wise
% likelihoods.
nSamples = length(p.Results.Angles);
nMixtures = length(p.Results.CProp);

compLik = zeros(nSamples, nMixtures);

% Compute component-wise likelihoods.
for idx = 1 : nMixtures
  compLik(:, idx) = (1 ./ ...
    (2 * pi * besseli(0, p.Results.Kappa(idx)))) .* ...
    exp(p.Results.Kappa(idx) .* cos(p.Results.Angles - p.Results.Mu(idx)));
end

% Accumulate component likelihoods to get the total likelihood.
likelihood = compLik * p.Results.CProp(:);

end

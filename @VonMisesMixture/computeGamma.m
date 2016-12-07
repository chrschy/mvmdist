function gamma = computeGamma(obj, angles, cProp, mu, kappa)
% COMPUTEGAMMA This function computes the responsibility matrix for a given
%   set of parameter values.
%
% REQUIRED INPUTS:
%   angles    - Nx1 vector, containing N angular values, ranged between -pi 
%               and pi.
%   cProp     - Kx1 vector, representing the proportions of all mixture 
%               components in the distribution.
%   mu        - Kx1 vector, containing the circular means of all K
%               components of the mixture distribution. The range of
%               each element is bounded in [-pi, pi].
%   kappa     - Kx1 vector, containing the concentration parameters
%               of the mixture distribution.  All elements are
%               nonnegative.
%
% OUTPUTS:
%   gamma     - NxK responsibility matrix.
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

% Get component-wise likelihoods
[~, cLik] = obj.computeLikelihood(p.Results.Angles, ...
  p.Results.CProp, p.Results.Mu, p.Results.Kappa);

% Compute gamma-matrix. Add a small positive number to gamma
% matrix to avoid numerical instabilities.
cGammas = bsxfun(@times, cLik, cProp') + sqrt(eps);

% This is to ensure that all elements of the gamma-matrix are
% greater than zero. Otherwise, problems can occur during the
% optimization.
gamma = bsxfun(@times, cGammas, 1 ./ sum(cGammas, 2));

end

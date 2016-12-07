function [mu, kappa, cProp] = estimateParameters(obj, angles, gamma)
% ESTIMATEPARAMETERS This function computes the maximum
%   likelihood estimates of the distribution parameters mu,
%   kappa and the mixing proportions for a given set
%   of angular values. The concentration parameter is estimated
%   using the approximation scheme introduced in [1].
%
% REQUIRED INPUTS:
%   angles - Nx1 vector, containing N angular values, ranged
%       between -pi and pi.
%   gamma - NxK responsibility matrix.
%
% OUTPUTS:
%   mu - Kx1 vector, containing the circular means of all K
%       components of the mixture distribution. The range of
%       each element is bounded in [-pi, pi].
%   kappa - Kx1 vector, containing the concentration parameters
%       of the mixture distribution.  All elements are
%       nonnegative.
%   cProp - Kx1 vector, representing the proportions of all
%       mixture components in the distribution.
%
% LITERATURE:
%   [1] D.J. Best, N.I. Fisher (1981): "The bias of the maximum
%       likelihood estimators of the von Mises-Fisher
%       concentration parameters"
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

p.addRequired('Gamma', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', '2d', 'nrows', length(angles)}) ...
  );

p.parse(angles, gamma);

% Compute mixing proportions.
cProp = sum(gamma) ./ size(gamma, 1);

% Ensure that each mixing proportion coefficient is greater than zero.
cProp = cProp + eps;
cProp = cProp ./ sum(cProp);

% Compute estimation parameters.
x = (gamma' * cos(angles)) ./ sum(gamma)';
y = (gamma' * sin(angles)) ./ sum(gamma)';

% Compute mean angles.
mu = atan2(y, x);

% Compute average distances.
d = sqrt(x.^2 + y.^2);

% Use approximation scheme from [1] to estimate kappa.
kappa = zeros(obj.nComponents, 1);

if any(d < 0.53)
  idx = d < 0.53;
  kappa(idx) = 2 .* d(idx) + d(idx).^3 + (5 .* d(idx).^5) ./ 6;
end

if any(0.53 <= d) && any(d < 0.85)
  idx = (0.53 <= d) & (d < 0.85);
  kappa(idx) = -0.4 + 1.39 .* d(idx) + 0.43 ./ (1 - d(idx));
end

if any(d >= 0.85)
  idx = d >= 0.85;
  kappa(idx) = 1 ./ (d(idx).^3 - (4 .* d(idx).^2) + (3 .* d(idx)));
end

% Force results to be column vectors.
cProp = cProp(:);

end

function likelihood = pdf(obj, angles)
% PDF This function returns the likelihood values for a set of angles.
%
% REQUIRED INPUTS:
%   angles          - Nx1 vector, containing N angular values, ranged
%                     between -pi and pi.
%
% OUTPUTS:
%   likelihood      - Nx1 vector containing the likelihood values for all
%                     input angles.
%   compLik         - NxK matrix, containing the component-wise likelihood
%                     values for all input angles.
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

p.parse(angles);

likelihood = obj.computeLikelihood(p.Results.Angles, ...
  obj.componentProportion, obj.mu, obj.kappa);

end

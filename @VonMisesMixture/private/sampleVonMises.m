function angles = sampleVonMises(mu, kappa, varargin)
% SAMPLEVONMISES This function generates samples from a unimodal von Mises
%   distribution. The sampling scheme used by this function was adopted 
%   from [1].
%
% REQUIRED INPUTS:
%   mu - Circular mean of the distribution. The parameter mu must be a real
%       valued scalar between -pi and pi.
%   kappa - Nonnegative, real-valued scalar concentration parameter of the
%       distribution.
%
% OPTIONAL INPUTS:
%   nSamples - Number of samples that should be generated. (default = 1)
%
% OUTPUTS:
%   angles - nSamples x 1 vector containing the sampled angular values.
%
% REFRENCES:
%   [1] L. Barabesi (1995): "Generating von Mises Variates by the
%       Ratio-of-Uniforms Method"
%
% AUTHOR:
%   Copyright (c) 2015      Christopher Schymura
%                           Cognitive Signal Processing Group
%                           Ruhr-Universitaet Bochum
%                           Universitaetsstr. 150
%                           44801 Bochum, Germany
%                           E-Mail: christopher.schymura@rub.de
%
% LICENSE INFORMATION:
%   This program is free software: you can redistribute it and/or
%   modify it under the terms of the GNU General Public License as
%   published by the Free Software Foundation, either version 3 of the
%   License, or (at your option) any later version.
%
%   This material is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <http://www.gnu.org/licenses/>

% Check inputs
p = inputParser();
defaultNSamples = 1;

p.addRequired( 'mu', @(x) validateattributes(x, {'numeric'}, ...
    {'real', 'scalar', '>=', -pi, '<=', pi}) );
p.addRequired( 'kappa', @(x) validateattributes(x, {'numeric'}, ...
    {'real', 'scalar', 'nonnegative'}) );
p.addOptional( 'nSamples', defaultNSamples, @(x) validateattributes(x, ...
    {'numeric'}, {'real', 'scalar', 'nonnegative'}) );
p.parse( mu, kappa, varargin{:} );

% Initialize sampling parameter
if p.Results.kappa > 1.3
    samplingParameter = 1 / sqrt( p.Results.kappa );
else
    samplingParameter = pi * exp( -p.Results.kappa );
end

% Allocate output
angles = zeros( p.Results.nSamples, 1 );

% Perform sampling
for idx = 1 : p.Results.nSamples
    while ( true )
        % Generate two random samples from the uniform
        % distribution
        randomSamples = rand( 2, 1 );
        
        % Compute angular value
        angle = ( samplingParameter * ...
            (2 * randomSamples(2) - 1) ) / randomSamples( 1  );
        
        % Check if angle is in valid range.
        if abs( angle ) > pi
            continue;
        end
        
        % Check first stopping condition
        if p.Results.kappa * angle^2 < 4 - 4 * randomSamples( 1 )
            break;
        end
        
        % Check second stopping condition
        if p.Results.kappa * cos( angle ) < ...
                2 * log( randomSamples(1) ) + p.Results.kappa
            continue;
        else
            break;
        end
    end
    
    % Assign sample value to output
    angles( idx ) = atan2( sin(angle + p.Results.mu), ...
        cos(angle + p.Results.mu) );
end

end
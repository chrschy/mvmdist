function angles = random(obj, varargin)
% RANDOM This function generates samples from the specified
%   mixture of von Mises distributions.
%
% OPTIONAL INPUTS:
%   nSamples        - Number of samples to be generated.
%   useMex          - Boolean variable, indicating if the mex-
%                     version of random number sampling should
%                     be used.
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

defaultNSamples = 1;
defaultUseMex = true;

p.addOptional('NSamples', ...
  defaultNSamples, ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'scalar', 'nonnegative'}) ...
  );

p.addOptional('UseMex', ...
  defaultUseMex, ...
  @islogical );

p.parse(varargin{:});

% Initialize resulting vector of angular samples.
angles = zeros(p.Results.NSamples, 1);

if p.Results.UseMex
  % Generate vector of component indices for sampling according to the
  % component proportion assigned to the model.
  compIdxs = 1 : obj.nComponents;
  compIdxVec = compIdxs(sum(bsxfun(@ge, rand(p.Results.NSamples, 1), ...
    cumsum(obj.componentProportion ./ sum(obj.componentProportion))'), 2) + 1 );
  
  currentSampleIdx = 1;
  for idx = 1 : obj.nComponents
    % Get number of samples for current component and sample using the mex
    % file for random number generation.
    nSamplesComp = sum(compIdxVec == idx);
    
    % Sample from current component.
    if nSamplesComp ~= 0      
      angles(currentSampleIdx : currentSampleIdx + nSamplesComp - 1) = ...
        sampleVonMisesMex(obj.mu(idx), obj.kappa(idx), nSamplesComp);
      
      currentSampleIdx = currentSampleIdx + nSamplesComp;
    end
  end
else
  % A random component index is picked during sampling, based on the
  % probabilities given by the individual weights. Therefore, the 
  % cumulative distribution has to be computed in advance.
  cumProbs = cumsum(obj.componentProportion);
  
  for idx = 1 : p.Results.NSamples
    % Select component to sample from and perform sampling from the 
    % selected component.
    cIdx = 1 + sum(cumProbs < rand());
    
    angles(idx) = sampleVonMises(obj.mu(cIdx), ...
      obj.kappa(cIdx));
  end
end

end

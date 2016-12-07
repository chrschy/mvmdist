classdef VonMisesMixture
  % VONMISESMIXTURE This class can be used to model a mixture of von
  %   Mises distributions. It provides functions to estimate the
  %   distribution parameters, compute likelihoods for given observations
  %   and draw samples from the distribution. Furthermore, it can be used
  %   for clustering applications involving circular data. Please refer
  %   to the individual function headers for further details.
  %
  % AUTHOR:
  %   Copyright (c) 2016      Christopher Schymura
  %                           Cognitive Signal Processing Group
  %                           Ruhr-Universitaet Bochum
  %                           Universitaetsstr. 150
  %                           44801 Bochum, Germany
  %                           E-Mail: christopher.schymura@rub.de
  
  properties (GetAccess = public, SetAccess = protected)
    % Number of mixture components.
    nComponents
    
    % Kx1 vector, containing the circular means of all k components of the 
    % mixture distribution. The range of each element is bounded in 
    % [-pi, pi].
    mu
    
    % Kx1 vector, containing the concentration parameters of the mixture 
    % distribution. All elements are nonnegative.
    kappa                   
    
    % Kx1 vector, representing the proportions of all mixture components in
    % the distribution.
    componentProportion
    
    % Log-likelihood achieved during parameter estimation.
    logLikelihood
    
    % Number of EM iterations performed for parameter estimation.
    nIterations
    
    % This flag indicates if the EM algorithm has converged to a local 
    % optimum during parameter estimation.
    hasConverged            
  end
  
  methods (Static, Hidden)
    [idx, centers] = ckmeans(angles, nClusters, varargin)
    
    [likelihood, compLik] = computeLikelihood(angles, weights, ...
      cProp, mu, kappa)
  end
  
  methods (Access = public)
    function obj = VonMisesMixture(varargin)
      % VONMISESMIXTURE This class constructor can be used to either 
      %   instantiate a "blank" distribution without any assigned 
      %   parameters or use additional input arguments to initialize with a
      %   given set of parameters.
      %
      % OPTIONAL INPUTS:
      %   componentProportion     - Kx1 vector representing the
      %                             proportions of all mixture
      %                             components in the distribution.
      %   mu                      - Kx1 vector, containing the
      %                             circular means of all K
      %                             components of the mixture
      %                             distribution. The range of each
      %                             element must be bounded in
      %                             [-pi, pi].
      %   kappa                   - Kx1 vector, containing the
      %                             concentration parameters of the
      %                             mixture distribution. All
      %                             elements must be nonnegative.
      
      p = inputParser();
      defaultComponentProportion = 1;
      defaultMu = 0;
      defaultKappa = 1;
      
      p.addOptional('ComponentProportion', ...
        defaultComponentProportion, ...
        @(x) validateattributes(x, ...
        {'numeric'}, ...
        {'real', 'vector'}) ...
        );
      
      p.addOptional('Mu', ...
        defaultMu, ...
        @(x) validateattributes(x, ...
        {'numeric'}, ...
        {'real', 'vector','>=', -pi, '<=', pi, ...
        'numel', length(varargin{1})}) ...
        );
      
      p.addOptional('Kappa', ...
        defaultKappa, ...
        @(x) validateattributes(x, ...
        {'numeric'}, ...
        {'real', 'vector', 'nonnegative', ...
        'numel', length(varargin{1})}) ...
        );
      
      p.parse(varargin{:});
      
      % Check if component proportions sum to one.
      if abs(sum(p.Results.ComponentProportion) - 1) > sqrt(eps)
        error('Component proportions must sum to one.');
      end
      
      obj.componentProportion = p.Results.ComponentProportion(:);
      obj.mu = p.Results.Mu(:);
      obj.kappa = p.Results.Kappa(:);
      
      % Get number of mixture components from specified component
      % proportions.
      obj.nComponents = length(p.Results.ComponentProportion);
    end
    
    likelihood = pdf(obj, angles)
    
    angles = random(obj, varargin)
    
    obj = fit(obj, angles, nComponents, varargin)
    
    [idx, nLogLik, posteriors] = cluster(obj, angles)
  end
  
  methods (Access = private)
    gamma = computeGamma(obj, angles, weights, cProp, mu, kappa)
    
    [mu, kappa, cProp] = estimateParameters(obj, angles, weights, gamma)
  end
end
%% SampleParamSets.m:

% FUNCTION NAME:
%  SampleParamSets
%
% DESCRIPTION:
%   Draw priors from parameter ranges.
%
% INPUTS:
%   nDraws: Integer, number of draws from the ranges.
%   parameter_priors: array of parameter ranges.
%
% OUTPUT:
%   Arrays containing the parameter sets.
%

function ParamSets = SampleParamSets(nDraws,parameter_priors)

nParams = size(parameter_priors, 1);

ParamSets = zeros(nParams,nDraws);
for i = 1:nParams
    
    minP = parameter_priors(i,1);
    maxP = parameter_priors(i,2);
    rangeP = maxP-minP;

    ParamSets(i,:) = ...
        minP + rangeP*lhsdesign(1,nDraws,'criterion','correlation');
    
end

end
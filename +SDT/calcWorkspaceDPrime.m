% [datROC,dP,optMet] = calcWorkspaceDPrime(TD,dataType)
%
% Collects the optimization algorthim data for the projection 2D workspace 
% or neural data and then calculates the d' value for the data.
%
% Inputs: 
%   TD          Trajectory Data structure 
%   dataType    Sets which type of data we are calculating d' (cursor or
%                   neural space)
%
% Outputs: 
%   dP        A single scalar value of d' 
%   datROC    A structure that returns the data necessary to remake an
%                ROC curve for the data, ie P, N, sensitivity (TPR = TP/P),
%                specificity (TN/N = 1-FPR), criterion, and t-test results
%                comparing the two distributions.
%   optMet    Structure with parameterized data passed to the optimization.
%
% Author:       Erinn Grigsby
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

function [datROC,dP,optMet] = calcWorkspaceDPrime(TD,dataType)
[~,~,condition] = unique([TD.targPos]','rows');

data = cell(1,length(TD));
idx = ones(length(TD),1);
for n = 1:length(TD)
    % Determine which space the d' should be calculated in
    if ismember(dataType,'pos') % The workspace
        data{n} = TD(n).brainKin.pos(:,1:2)';
    elseif ismember(dataType,'xorth') % The orthonormalized neural space
        data{n} = TD(n).GPFA.xorth;
    elseif ismember(dataType,'xsm') % The neural space
        data{n} = TD(n).GPFA.xsm;
    end

    if isempty(data{n})
        idx(n) = 0;
    end
end

% Remove the empty index
idx = idx & [TD.successful]';
data = data(idx==1);
condition = condition(idx==1);

[optMet,~] = opt.get_data(data,condition);

% Create the midpoint normalization vector
vecStrt = optMet.mu_AB - optMet.mu_BA;
vecStrt = vecStrt./norm(vecStrt);

% Project data on to the vector
midPt = optMet.x_mid;
parMidPt = vecStrt'*midPt;
a2D = parMidPt(condition==1);
b2D = parMidPt(condition==2);

% Calculate the d' for the projection of the midpoint
[dP,datROC] = SDT.calcSensitivityIdx(a2D',b2D');
end
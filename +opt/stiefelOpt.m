function M = stiefelOpt(P,objFn,gradFn,varargin)
% Perform gradient descent over the Stiefel manifold
%
% Inputs:
%   P       data
%   objFn   Handle of objective function to evaluate
%   gradFn  Handle of gradient function to evaluate
%
% Provide input data in a structure.  This way it can be passed to the
% objective and gradient functions in a generic form
%
% Optional Inputs:
%   plotFigs        Plot convergence of objective function
%   verbose         Display objective function values when fitting
%   objFnParams     Objective function parameters (weights of the various terms)
%   R               Number of orthonormal projection vectors to return
%                       (currently hard-coded)
%
% Outputs:
%   M       Projection matrix
%
% Implementation based on Cunningham and Ghahramani 2015.
%
% Created by Alan Degenhart
% Copyright (C) by Erinn Grigsby and Alan Degenhart
% Emails: erinn.grigsby@gmail.com or alan.degenhart@gmail.com

% Parse optional arguments
plotFigs = false;
verbose = true;     % Display objective function values when fitting
objFnParams = [];   % Objective function parameters (weights of the various terms)
R = 2;              % Number of orthonormal projection vectors to return (currently hard-coded)
assignopts(who,varargin);

% Parameters
maxIter = 1000;     % Maximum number of iterations
D = P.D;            % Number of dimensions in full space
nRestarts = 5;      % Number of random restarts

deltaB = 0.9;
epsF = 1e-10;       % Convergence criteria

O = repmat(struct(),nRestarts,1);
for i = 1:nRestarts
    % Initialize M (randomly)
    A = randn(D,R);
    M = orth(A);    % Columns of M will be orthonormal

    if verbose
        fprintf('\n\nPerforming optimization ... \n')
    end

    % Perform gradient descent iterations
    J = nan(maxIter+1,1);   % '+1' to account for the initial random value of M
    [J(1),JTtemp,objStr] = objFn(M,P,objFnParams);
    JT = nan(maxIter+1,length(JTtemp));
    JT(1,:) = JTtemp;
    converged = false;
    currIter = 1;
    while (~converged && currIter <= maxIter)
        Mprev = M;
        b = 0.1; % Step size (beta)

        % Step 1: Calculate free gradient
        Z = gradFn(M,P,objFnParams);

        % Step 2: Compute search direction
        if R>2
            Z(:,3:R) = 0;
        end
        Z = searchDir(-Z,M);
        
        % Step 3: Line search with retraction
        fM = objFn(M,P,objFnParams);
        fR = objFn(retract(b*Z,M),P,objFnParams);
        dF = fM - fR;
        maxLsIter = 500;
        dfMat = nan(maxLsIter,1);
        lsIter = 1;
        while dF < 0
            % Adjust B
            b = b * deltaB;

            % Evaluate step
            fM = objFn(M,P,objFnParams);
            fR = objFn(retract(b*Z,M),P,objFnParams);
            dF = fM - fR;
            dfMat(lsIter) = dF;

            % Debug to check for line search convergence
            if lsIter == maxLsIter
                keyboard
            end
            lsIter = lsIter + 1;
        end

        % Step 4: Update estimate of M
        M = retract(b*Z,M);

        % Step 5: Determine convergence (difference between gradient is less
        % than some convergence criteria)
        currIter = currIter + 1;
        [J(currIter),JT(currIter,:)] = objFn(M,P,objFnParams);
        dJ = J(currIter) - J(currIter - 1);
        converged = abs(dJ) < epsF;

        % Display convergence info
        if (mod(currIter,1) == 0) && verbose
            fprintf('Iter: %d, J = %0.5g, dJ = %0.5e\n',currIter,J(currIter),dJ)
        end
    end
    
    % Put results from current iteration into 'O' structure
    O(i).M = M;
    O(i).J = J;
    O(i).dJ = dJ;
    O(i).JT = JT;
    O(i).Jconv = J(currIter);
    O(i).converged = converged;
    O(i).iter = currIter;
end

% Select best projection across random restarts
allJconv = [O.Jconv];
[~,I] = min(allJconv);
M = O(I).M;
J = O(I).J;
JT = O(I).JT;

% Plot figures if desired
if plotFigs
    % Plot convergence of objective function (and individual terms)
    plot(J,'k'); hold on
    legStr{1} = 'obj Fn';
    for i = 1:size(JT,2)
        plot(JT(:,i),'LineStyle','--')
    end
    legStr = {'Obj Fn.'};
    legStr = [legStr objStr];
    L = legend(legStr);
    set(L,'Interpreter','latex')
    
    % Evaluate objective function for original projection
    M_opt = zeros(D,R);
    M_opt(1,1) = 1;
    M_opt(2,2) = 1;
    J_opt = objFn(M_opt,P,objFnParams);
    set(gca,'box','off','TickDir','out')
    xlabel('Iteration')
    ylabel('Objective Fn. Value')
    title('Optimization Convergence')

    % Project P into identified projection and plot
    P_opt = P;
    P_opt.mu_A = M' * P.mu_A;
    P_opt.mu_B = M' * P.mu_B;
    P_opt.mu_AB = M' * P.mu_AB;
    P_opt.mu_BA = M' * P.mu_BA;
    P_opt.sig_AB = M' * P.sig_AB * M;
    P_opt.sig_BA = M' * P.sig_BA * M;
end

end % EOF

% Search direction
function Z = searchDir(Z,M)
    D = size(M,1);
    SK = (1/2) * (M'*Z - Z'*M);
    Z = M * SK + (eye(D) - M*M') * Z;
end

% Retraction
function Z = retract(Z,M)
    R = size(M,2);
    Z = (M + Z)*(eye(R) + Z'*Z)^(-1/2);
end
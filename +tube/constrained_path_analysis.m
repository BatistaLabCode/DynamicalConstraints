% el.tube.constrained_path_analysis  Perform analysis for constrained path
% experiment.


function [IT, T] = constrained_path_analysis(IT, varargin)

% Optional agruments
intTargNum = 1;
blockSize = 50;
kinSource = 'brain';
tubeSegment = 0;

assignopts(who,varargin);

% Check for valid washout data.  If this exists, tube constraints will be
% applied to these trials as well.
validWashout = false;
if ~isempty(IT.TDunconstrainedWashout)
    validWashout = true;
end

% Get number of tube constraints to apply.
nTube = length(IT.TDconstrained);

% Find number of start positions.  Here we use the first constrained tube
% block, as sometimes only one start position is 
TD = IT.TDconstrained{1};
startPos = [TD.startPos]';
uniStartPos = unique(startPos,'rows');
nStartPos = size(uniStartPos,1);

% Get unique start positions and associated tubes
TUall = repmat(TubeObject(),nTube,nStartPos);
tubeTrialsAll = nan(nTube,nStartPos);
tubeTrialsSuc = nan(nTube,nStartPos);
tubeRadius = nan(nTube,nStartPos);
for i = 1:nTube
    % Get trajectory and tube data, set tube code.  Currently tube
    TDtemp = IT.TDconstrained{i};
    TU = [TDtemp.tube];
    TU = TU.setTubeCode;
    
    % Find start positions and get tube objects for each position
    tubeStartPos = [TDtemp.startPos]';
    for j = 1:nStartPos
        % Get tube object
        startPosMask = ismember(tubeStartPos,uniStartPos(j,:),'rows');
        TUtemp = TU(startPosMask);
        TUall(i,j) = TUtemp(1);
        
        % Get tube radius
        r = TUtemp(1).window(:,1);
        r = r(~isnan(r));
        tubeRadius(i,j) = r(1);
        TUall(i,j).radius = r(1);
        
        % Calculate acutal success rate
        TDtube = TDtemp(startPosMask);
        tubeTrialsAll(i,j) = length(TDtube);
        tubeTrialsSuc(i,j) = sum([TDtube.successful]);
    end
end
tubeSR = tubeTrialsSuc./tubeTrialsAll;

% Apply tube constraints to unconstrained trials
[TD_unconst_tube,condTrialsAll,condTrialsSuc] = apply_tube_constraint( ...
    IT.TDunconstrained,uniStartPos,TUall,kinSource,tubeSegment);
condSR = condTrialsSuc./condTrialsAll;

% Apply tube constraints to unconstrained trials (washout)
if validWashout
    [TD_unconst_tube_wo,condTrialsAll_wo,condTrialsSuc_wo] = apply_tube_constraint( ...
        IT.TDunconstrainedWashout,uniStartPos,TUall,kinSource,tubeSegment);
    condSR_wo = condTrialsSuc_wo./condTrialsAll_wo;
else
    condSR_wo = [];
    TD_unconst_tube_wo = [];
end

% Add data to IT structure
IT.successRate = tubeSR;
IT.expectedSuccessRate = condSR;
IT.expectedSuccessRate_washout = condSR_wo;
IT.constrainedTubeRadius = tubeRadius;
IT.startTargPos = unique(startPos,'rows');

% Loop over tubes
tubeBlockSR = nan(nTube,2); % Early and late success rate
tubeBestSR = nan(nTube,1);  % Best success rate
TD_tube_best = cell(1,nTube);
for i = 1:nTube
    
    TD = IT.TDconstrained{i};
    
    % Calculate running success rate to find the best block of trials
    scMask = logical([TD.successful]);
    sR = conv(scMask,ones(1,blockSize),'full')/blockSize;
    sR = sR(1:length(scMask));
    if length(scMask) < blockSize
        tubeBestSR(i) = sum(scMask)/length(scMask);
        tubeBlockSR(i, 1) = tubeBestSR(i);
        I = length(scMask);
    else
        [tubeBestSR(i), I] = max(sR);
        tubeBlockSR(i, 1) = sum(scMask(1:blockSize))/blockSize;
    end
    
    % Get best block of trials
    if I < blockSize
        idx_best = 1:I;
    else
        idx_best = (I - blockSize + 1):I;
    end
    TD_tube_best{i} = TD(idx_best);  % Keep failed trials for now
end

%--------------------------------------------------------------------------
% Calculate observed and expected success rates

% Split unconstrained trials up into predicted difficult and predicted easy
% starting positions.
TD_unconst = IT.TDunconstrained;
start_pos = [TD_unconst.startPos]';
mask_difficult = ismember(start_pos, TD_tube_best{end}(1).startPos', 'rows');
TD_unconst_difficult = TD_unconst(mask_difficult);
TD_unconst_easy = TD_unconst(~mask_difficult);

% Calculate expected success rate for best 50 trials
r = 1:200;
sr_const_best = tube.calcExpectedSR(TD_tube_best{end},TUall(1),r,...
    'kinSource',kinSource);
sr_unconst = tube.calcExpectedSR(TD_unconst_difficult,TUall(1),r, ...
    'startPos', TD_tube_best{end}(1).startPos','kinSource',kinSource);

% Apply smallest tube to unconstrained trials
[~, sr_unconst_tube_mask] = tube.calcExpectedSR( ...
    TD_unconst_difficult, TUall(1), TUall(end).radius,'kinSource',kinSource);
TD_unconst_difficult_best = TD_unconst_difficult(sr_unconst_tube_mask{1});
if tubeSegment & ~isempty(TD_unconst_difficult_best)
    sr_unconst_best = tube.calcExpectedSR(TD_unconst_difficult_best,...
        TUall(1),r,'startPos', TD_tube_best{end}(1).startPos',...
        'kinSource',kinSource,'plotSegment',TUall(1).segment);
elseif ~isempty(TD_unconst_difficult_best)
    sr_unconst_best = tube.calcExpectedSR(TD_unconst_difficult_best,...
        TUall(1),r,'startPos', TD_tube_best{end}(1).startPos',...
        'kinSource',kinSource);
else
    sr_unconst_best = zeros(1,length(r));
end

if ~isempty(TD_unconst_easy)
    unconst_easy_flag = true;
    
    % Create a tube for the predicted easy direction, based on data
    % provided
    tubePath = nan(100,3);
    tubeStartPos = TD_unconst_easy(1).startPos';
    tubeStartPos(2,:) = TD_unconst_easy(1).targPos';
    for k = 1:3
        tubePath(:,k) = spline([1 100],tubeStartPos(1:2,k),1:100);
    end
    TU = TubeObject();
    TU.path = tubePath;
    TU.radius = r(1);
    TU.window = r(1);
                
    sr_unconst_easy = tube.calcExpectedSR(TD_unconst_easy, ...
        TU, r,'startPos', TD_unconst_easy(1).startPos','kinSource',kinSource);
else
    unconst_easy_flag = false;
    sr_unconst_easy = nan(1, length(r));
end

% Fit sigmoid to success rate and find tube radius where performance is 50%
[P_const] = fit_sigmoid(r, sr_const_best);
[P_unconst] = fit_sigmoid(r, sr_unconst);
[P_unconst_best] = fit_sigmoid(r, sr_unconst_best);

if unconst_easy_flag
    [P_unconst_easy] = fit_sigmoid(r, sr_unconst_easy);
else
    P_unconst_easy.x50 = nan;
    P_unconst_easy.x_fit = nan;
    P_unconst_easy.y_fit = nan;
end

% Add trajectory data to IT object
IT.intermediate_target_num = intTargNum;
IT.block_size = blockSize;
IT.TD_tube_best = TD_tube_best;
IT.TD_unconst_tube = TD_unconst_tube;
IT.TD_unconst_tube_washout = TD_unconst_tube_wo;
IT.tubeObject = TUall;

% General data
T.subject = IT.subject;
T.date = IT.date;
T.int_targ_num = intTargNum;
T.start_pos = uniStartPos;
T.tubeRadius = tubeRadius;
T.actualSR = tubeSR;
T.actualSR_best = tubeBestSR;
T.actualSR_block = tubeBlockSR;
T.expectedSR = condSR;

% Data used for sigmoid fits
T.r = r;
T.sr_unconst_all = sr_unconst;
T.sr_unconst_best = sr_unconst_best;
T.sr_const_best = sr_const_best;
T.sr_unconst_easy = sr_unconst_easy;

% Results of sigmoid fits
T.P_unconst = P_unconst;
T.P_unconst_best = P_unconst_best;
T.P_const = P_const;
T.P_unconst_easy = P_unconst_easy;

% Other info
T.int_targ_dist = norm(TD_tube_best{end}(1).targPos);  % Distance of intermediate target from the center

end % EOF

%--------------------------------------------------------------------------
function [TDtubeAll,condTrialsAll,condTrialsSuc] = ...
    apply_tube_constraint(TD,uniStartPos,TUall,kinSource,tubeSegment)
% Inputs
% - nTube
% - nStartPos
% - uniStartPos
% - TUall

nTube = size(TUall,1);
nStartPos = size(TUall,2);
srcStr = [kinSource 'Kin'];
plotSegment = [];

% Loop over tube radii and apply tubes to unconstrained trials
TDtubeAll = cell(nTube,nStartPos);
condTrialsAll = nan(nTube,nStartPos);
condTrialsSuc = nan(nTube,nStartPos);
for i = 1:nTube
    TDtube = TD;
    startPos = [TDtube.startPos]';
    validMask = false(1,length(TDtube));
    
    % Loop over start positions
    for j = 1:nStartPos
        % Get trials for start target and associated tube
        trialMask = ismember(startPos,uniStartPos(j,:),'rows');
        TDtemp = TDtube(trialMask);
        TU = TUall(i,j);
        
        if tubeSegment
            plotSegment = TU.segment;
        end
        
        % Loop over trials and apply tube constraint
        for k = 1:sum(trialMask)
            pos = TDtemp(k).(srcStr).pos;
            contMask = TU.calcContainment(pos,'plotSegment',plotSegment);
            
            % If not all points are contained in the tube, update the
            % trajectory data object
            if sum(contMask) ~= length(contMask)
                % Find index of first failed point
                offsetInd = find(~contMask,1,'first') - 1;
                
                % Truncate kinematics
                TDtemp(k).kinTime = TDtemp(k).kinTime(1:offsetInd);
                TDtemp(k).pos = TDtemp(k).pos(1:offsetInd,:);
                TDtemp(k).vel = TDtemp(k).vel(1:offsetInd,:);
                TDtemp(k).acc = TDtemp(k).acc(1:offsetInd,:);
                TDtemp(k).(srcStr).time = TDtemp(k).(srcStr).time(1:offsetInd);
                TDtemp(k).(srcStr).pos = TDtemp(k).(srcStr).pos(1:offsetInd,:);
                TDtemp(k).(srcStr).vel = TDtemp(k).(srcStr).vel(1:offsetInd,:);
                TDtemp(k).(srcStr).acc = TDtemp(k).(srcStr).acc(1:offsetInd,:);
                
                % Update success code
                TDtemp(k).successful = false;
            end
            
            % Update tube object
            TDtemp(k).tube = TU;
        end
        
        % Calculate success rate for condition
        condTrialsAll(i,j) = sum(trialMask);
        condTrialsSuc(i,j) = sum([TDtemp.successful]);
        
        % Update trajectory data
        validMask(trialMask) = true;
        TDtube(trialMask) = TDtemp;
    end
    % Get rid of any start positions that were not analyzed
    TDtube = TDtube(validMask);
    TDtubeAll{i} = TDtube;
end

end

%--------------------------------------------------------------------------
% Fit sigmoid to data.  This is used to analyze the success rate for 
function [P] = fit_sigmoid(x,y)
    % Define function to fit
    sigfunc = @(A, x)(A(1) ./ (1 + exp(-A(2)*x + A(3))));
    
    % Set initial conditions and fit
    A0 = ones(1,3); % Initial values
    A0(1) = max(y) - min(y);
    A0(3) = (max(x) - min(x))/2;
    A_fit = nlinfit(x, y, sigfunc, A0);
    
    % Predict values
    y_fit = sigfunc(A_fit, x);
    x_fit = x;
    
    % Calculate 50% point
    x50 = (log(2*A_fit(1) - 1) - A_fit(3))/(-A_fit(2));

    P.x_fit = x_fit;
    P.y_fit = y_fit;
    P.A = A_fit;
    P.x50 = x50;
end

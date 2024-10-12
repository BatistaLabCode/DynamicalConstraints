%% Load the data and the decoder information, and removed the unnecessary
%
function [invalid] = batch_analysis_rand_Projection(D,varargin)
N = 250; % Number of projection iterations
n_perm = 5; % Number of random perturbations
savePath = [];
n_cores = 4;
tasks = {'tt_int','tt_rot'};
assignopts(who, varargin);

% Determine where to save data
if isempty(savePath)
    savePath = uigetdir();
end

% Load dataset info
if isempty(D)
    [File_D, loc_D] = uigetfile('*.mat','Load the file containing the summary of sessions');
    load([loc_D File_D])
end

% Turn on the parallel pool process
parpool(n_cores)
dir_list = db.get_task_datasets(D,tasks);
invalid = []
for n = 1:size(dir_list,1)
    try
        tFF = calc_random_flowField(dir_list(n,1),dir_list(n,2),savePath,...
            'N',N,'n_perm',n_perm);
        fprintf('Session %03d of %03d completed \n',n,size(dir_list,1))
    catch
        sprintf('Issue with session %03d of %03d completed...skipping\n',n,size(dir_list,1))
        invalid = [invalid n]
    end
end
end

function [FF] = calc_random_flowField(intuitive_dir,rotated_dir,savePath,varargin)
N = 500; % Number of projection iterations
n_perm = 30; % Number of random perturbations

assignopts(who,varargin);

subject = intuitive_dir.subject;
session = intuitive_dir.dataset;

% Load the trajectory data in
intTD = TrajectoryData().load(intuitive_dir.trajectory{1});
rotTD = TrajectoryData().load(rotated_dir.trajectory{1});
uniTarg = unique([rotTD.targPos]','rows');
centerPos = mean(uniTarg);

intName = intuitive_dir.translated{1};
rotName = rotated_dir.translated{1};
intPos = regexp(intName,'_');
intPos = [intPos(1)+1:intPos(2)-1];
rotPos = regexp(rotName,'_');
rotPos = [rotPos(1)+1:rotPos(2)-1];

% Normalize the trajectory data
intTD = intTD.normalize(centerPos);
rotTD = rotTD.normalize(centerPos);

% load the intuitve decoder
P_int = load(fullfile(intuitive_dir.base,intTD(1).decoderName));
P_int = P_int.bci_params;

% Load the results
decodeNum = P_int.decoderID;
method = P_int.decoderType;
xDim = P_int.xDim;
runDir = sprintf('mat_results/run%03d', decodeNum);
fname = sprintf('%s/%s_xDim%02d', runDir, method, xDim);
result = load(fullfile(intuitive_dir.base,'analysis',fname));

% Create the individual session save folder for the iterations
savePath = fullfile(savePath,sprintf('Iter_%04d_Perm_%01d',N,n_perm));
mkdir(savePath)
savePathSess = fullfile(savePath,'sessions');%,[subject session]);
mkdir(savePathSess)

% Calculate the projections
nDim = 10;
nR = 2;
rand('state',0) % Set the random stat so that we always getting the same initial projections
M = randn(nDim,nR,N);

% TODO: Determine the minimum number of trials to use. Assume 75% of the
% smallest group presented.
numTrl = 20;

%% Calculate the permutations
uniTarg = unique([rotTD.targPos]','rows');
permD = repmat(struct('cond',[],'target',[],'idx1',[],'idx2',[]),...
    2*size(uniTarg,1),1);
for n = 1:size(uniTarg,1)
    maskI = ismember([intTD.targPos]',uniTarg(n,:),'rows');
    maskR = ismember([rotTD.targPos]',uniTarg(n,:),'rows');
    
    posI = find(maskI);
    posR = find(maskR);
    
    intIdx = util.getPermutationSamples(sum(maskI),'numIter',n_perm,'numTrl',numTrl);
    rotIdx = util.getPermutationSamples(sum(maskR),'numIter',n_perm,'numTrl',numTrl);
    
    % Add the intuitive data
    permD(n).cond = 'int';
    permD(n).target = round(uniTarg(n,:));
    permD(n).idx1 = posI([intIdx(:,1).idx]');
    permD(n).idx2 = posI([intIdx(:,2).idx]');
    
    % Add the rotated data
    permD(n+size(uniTarg,1)).cond = 'rot';
    permD(n+size(uniTarg,1)).target = round(uniTarg(n,:));
    permD(n+size(uniTarg,1)).idx1 = posR([rotIdx(:,1).idx]');
    permD(n+size(uniTarg,1)).idx2 = posR([rotIdx(:,2).idx]');
end

% Save the perturbation and projection data
saveName = sprintf('%s%s_int%s_rot%s_randProjectionInfo',...
    subject,session,intName(intPos),rotName(rotPos));
save(fullfile(savePathSess,saveName),'M','permD')

%% Calculate the GPFA space: Load the intuitive decoder and apply to the
% intuitive and rotated mapping sessions
[~,~,TT] = orthogonalize(zeros(P_int.xDim,1),result.estParams.C);
optInput = {'spikeCountSrc','decodeSpikeCounts','useT',1,'TT',TT,'trunGPFA',1};
[intTD] = predictDecodeState_GPFA(intTD,P_int,optInput{:});
[rotTD] = predictDecodeState_GPFA(rotTD,P_int,optInput{:});%'spikeCountSrc','decodeSpikeCounts','useT',1)

%% Iterate through and create the flow field. NOTE: This step uses parallel 
% pool processing (comment: 'parfor n = 1:N' and uncomment: for n = 1:N to
% turn off parallel processing)
iA = [permD(1:2).idx1];
iB = [permD(1:2).idx2];
rA = [permD(3:4).idx1];
rB = [permD(3:4).idx2];

parfor n = 1:N
%for n = 1:N
% Find the projection and apply it.
M1 = M(:,:,n);
M1 = orth(M1);
proj_int_TD = opt.applyProjection(intTD,M1,'savePos',1);
proj_int_rev_TD = reverseDirection(proj_int_TD);
proj_rot_TD = opt.applyProjection(rotTD,M1,'savePos',1);
proj_rot_rev_TD = reverseDirection(proj_rot_TD);

% Iterate throug the different groups
[ii1,ii2,iv1,iv2,ia1,ia2,rr1,rr2,rv1,rv2,ra1,ra2,ir1,ir2] = deal(...
    repmat(struct('condition',[],'n_valid',[],'angValid',[],'angMean',[],...
    'angMedian',[],'magValid',[],'magMean',[],'magMedian',[],'mseValid',...
    [],'mseMean',[],'mseMedian',[]),n_perm,1));

for k = 1:n_perm
    % Calculate intuitive maps A and B
    int_A = el.flow.calc_flow_field(proj_int_TD(iA(k,:)),1,'intuitive_A',...
        'min_pts_per_voxel',2,'ax_lim',[-20 20]);
    int_B = el.flow.calc_flow_field(proj_int_TD(iB(k,:)),1,'intuitive_B',...
        'min_pts_per_voxel',2,'ax_lim',[-20 20]);
    int_B_rev = el.flow.calc_flow_field(proj_int_rev_TD(iB(k,:)),1,...
        'reverse_intuitive_B','min_pts_per_voxel',2,'ax_lim',[-20 20]);
    
    % Calculate rotated maps A and B
    rot_A = el.flow.calc_flow_field(proj_rot_TD(rA(k,:)),1,'rotated_A',...
        'min_pts_per_voxel',2,'ax_lim',[-20 20]);
    rot_B = el.flow.calc_flow_field(proj_rot_TD(rB(k,:)),1,'rotated_B',...
        'min_pts_per_voxel',2,'ax_lim',[-20 20]);
    rot_B_rev = el.flow.calc_flow_field(proj_rot_rev_TD(rB(k,:)),1,...
        'reverse_rotated_B','min_pts_per_voxel',2,'ax_lim',[-20 20]);
    
    % Separate the trajectories by target pair
    fi_a = separateFFconditions(int_A);
    fi_b = separateFFconditions(int_B);
    fir_b = separateFFconditions(int_B_rev);
    fr_a = separateFFconditions(rot_A);
    fr_b = separateFFconditions(rot_B);
    frr_b = separateFFconditions(rot_B_rev);
    
    % Run the comparisions: 1-2) fi_a/fi_b[2] 3-4) fi_a/fir_b[2], 5-6)
    % fi_a/fi_a(alt) 7-8) fr_a/fr_b[2] 9-10) fr_a/frr_b[2] 11-12) fr_a/fr_a(alt)
    % 13-14) fi_a/fr_a[2]

    % 1-2)
    ii1(k) = flow.compare_flow_fields(fi_a(1), fi_b(1),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/int (1)',n,k));
    ii2(k) = flow.compare_flow_fields(fi_a(2), fi_b(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/int (2)',n,k));

    % 3-4)
    iv1(k) = flow.compare_flow_fields(fi_a(1), fir_b(1),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/int_rev  (1)',n,k));
    iv2(k) = flow.compare_flow_fields(fi_a(2), fir_b(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/int_rev (2)',n,k));

    % 5-6)
    ia1(k) = flow.compare_flow_fields(fi_a(1), fi_b(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/int_alt',n,k));
    ia2(k) = flow.compare_flow_fields(fi_b(1), fi_a(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/int_alt',n,k));
    
    % 7-8)
    rr1(k) = flow.compare_flow_fields(fr_a(1), fr_b(1),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, rot/rot (1)',n,k));
    rr2(k) = flow.compare_flow_fields(fr_a(2), fr_b(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, rot/rot (2)',n,k));

    % 9-10)
    rv1(k) = flow.compare_flow_fields(fr_a(1), frr_b(1),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, rot/rot_rev (1)',n,k));
    rv2(k) = flow.compare_flow_fields(fr_a(2), frr_b(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, rot/rot_rev (2)',n,k));

    % 11-12)
    ra1(k) = flow.compare_flow_fields(fr_a(1), fr_b(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, rot/rot_alt',n,k));
    ra2(k) = flow.compare_flow_fields(fr_b(1), fr_a(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, rot/rot_alt',n,k));
    
    % 13-14)
    ir1(k) = flow.compare_flow_fields(fi_a(1), fr_a(1),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/rot (1)',n,k));
    ir2(k) = flow.compare_flow_fields(fi_a(2), fr_a(2),'saveSimp',1,...
        'cond_str', sprintf('Iter: %d, Perm: %d, int/rot (2)',n,k));
end

FlowIteration = table(ii1,ii2,iv1,iv2,ia1,ia2,rr1,rr2,rv1,rv2,ra1,ra2,ir1,ir2);

% Save the iteration data - Optional if you want to save the individual
% iterations. Not necessary since it is possible to recreate the flow field
% from the other saved information.
% saveName = sprintf('%s%s_int%s_rot%s_flowRandProject_iter%03d',...
%     subject,session,intName(intPos),rotName(rotPos),n);
% save(fullfile(savePathSess,saveName),'FlowIteration','M1')

% Save the mean data
measField = {'angMean','angMedian','magMean','magMedian','mseMean','mseMedian','n_valid'};
condField = {'ii','iv','ia','ir','rr','rv','ra',};
for k = 1:length(measField)
    FF(n,1).(measField{k}).subject = subject;
    FF(n,1).(measField{k}).dataset = session;
    FF(n,1).(measField{k}).targPair = uniTarg;
    for m = 1:length(condField)
        FF(n,1).(measField{k}).(condField{m})(1) = nanmean([FlowIteration.([condField{m} '1']).(measField{k})]);
        FF(n,1).(measField{k}).(condField{m})(2) = nanmean([FlowIteration.([condField{m} '2']).(measField{k})]);
    end
end
end

% Save the summary flow field information.
saveName = sprintf('%s%s_int%s_rot%s_flowRandProject_results',...
    subject,session,intName(intPos),rotName(rotPos));
save(fullfile(savePath,saveName),'FF')
end

%% Functions
function [revTD] = reverseDirection(TD)
% Reverse the direction of the position data, need for calculating a
% portion of the flow field.
revTD = TD;

for n = 1: size(TD,1)
    revTD(n).pos = flipud(TD(n).pos);
end
end

function [FFsep] = separateFFconditions(FF)

nTarg = size(FF.kin.start_pos,1);

FFsep = repmat(FF,nTarg,1);

for n = 1:nTarg
    FFsep(n).condition = [FFsep(n).condition '_target' num2str(n)];
    
    % Simplify the kin data
    FFsep(n).kin.start_pos = FFsep(n).kin.start_pos(n);
    FFsep(n).kin.all_pos = FFsep(n).kin.all_pos(n);
    FFsep(n).kin.all_vel = FFsep(n).kin.all_vel(n);
    FFsep(n).kin.hist_counts = FFsep(n).kin.hist_counts(n);
    FFsep(n).kin.d_peak = FFsep(n).kin.d_peak(n);
    
    % Simplify the grid data
    FFsep(n).grid.num_pts_start_pos = FFsep(n).grid.num_pts_start_pos(:,:,n);
    FFsep(n).grid.num_pts = FFsep(n).grid.num_pts_start_pos;
    FFsep(n).grid.VX_start_pos = FFsep(n).grid.VX_start_pos(:,:,n);
    FFsep(n).grid.VX = FFsep(n).grid.VX_start_pos;
    FFsep(n).grid.VY_start_pos = FFsep(n).grid.VY_start_pos(:,:,n);
    FFsep(n).grid.VY = FFsep(n).grid.VY_start_pos;
    FFsep(n).grid.grid_cts = squeeze(FFsep(n).grid.grid_cts(1,:,:));
    FFsep(n).grid.grid_hist = FFsep(n).grid.grid_hist(1,:);
end

end
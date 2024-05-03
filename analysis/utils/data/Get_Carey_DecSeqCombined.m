function out = Get_Carey_DecSeqCombined(cfg_in)
% function out = Get_Carey_DecSeqCombined(cfg_in)
%
% run decoding sequence detection analysis from within session data folder
%
% assumes cfg_in.output_dir exists; image and data files are written into this
%
% CONFIGS:
%
% cfg_def = [];
% cfg_def.output_dir = 'files';
% cfg_def.output_file_prefix = []; % when writing files, prefix this
% cfg_def.dt = 0.025; % this is used for a firing rate threshold on Q matrix -- is that important?
% cfg_def.TCsmooth = 1; % bins SD; vary this and the next
% cfg_def.QsmoothSD = 0.002; % time SD; vary this with previous
% cfg_def.minSpikes = 25; % remove cells with less than this number of spikes during run
% cfg_def.nBins = 112; % ~3cm bins for full maze (334 cm); NOTE this is number of bin edges
% cfg_def.trialStartOffset = -1; % start "run" at this time relative to center pb break (in s)
% cfg_def.tc_baseline = 0.1; % baseline firing rate, replaces zeros in TC when unsmoothed with smoothed Q
% cfg_def.load_questionable_cells = 1;
% cfg_def.includeAllCells = 1; % otherwise, only place cells
% cfg_def.trackExcludeStart = 20; % exclude this amount (in cm) from start and end of track
% cfg_def.trackExcludeEnd = 15;
% cfg_def.nSpikesHist = -0.5:105;
% cfg_def.nMinNeurons = 4; % minimum number of neurons that needs to be active
% cfg_def.maxJump_cm = 40;
% cfg_def.minSeqLength = 10; % note, this is in bins
% cfg_def.nMaxNanSkipSequential = 0; % NaNs that can be skipped without breaking sequence
% cfg_def.plotOutput = 0;
% cfg_def.Qdt = 0.005; % this is the binsize used for decoding
% cfg_def.Qboxcar = 5; % boxcar smoothing of spike counts (in bins). So Qdt = 0.005 and Qboxcar = 5 give a true bin size of 25ms, moved in 5ms steps.
% cfg_def.writeFiles = 1;
% cfg_def.removeInterneurons = 0;
% cfg_def.postCPonly = 1; % if 1, only include maze arms (beyond choice point)
%
% MvdM 2015

%% master config
cfg_def = [];
cfg_def.output_dir = 'files';
cfg_def.output_file_prefix = []; % when writing files, prefix this
cfg_def.dt = 0.1; % either vary this (keeping others constant), or
cfg_def.TCsmooth = 1; % bins SD; vary this and the next
cfg_def.QsmoothSD = 0; % time SD; vary this with previous
cfg_def.minSpikes = 25; % remove cells with less than this number of spikes during run
cfg_def.nBins = 112; % ~3cm bins for full maze (334 cm); NOTE this is number of bin edges
cfg_def.trialStartOffset = -1; % start "run" at this time relative to center pb break (in s)
cfg_def.tc_baseline = 0.1; % baseline firing rate, replaces zeros in TC when unsmoothed with smoothed Q
cfg_def.load_questionable_cells = 1;
cfg_def.includeAllCells = 1; % otherwise, only place cells
cfg_def.trackExcludeStart = 20; % exclude this amount (in cm) from start and end of track
cfg_def.trackExcludeEnd = 15;
cfg_def.nSpikesHist = -0.5:105;
cfg_def.nMinNeurons = 4;
cfg_def.maxJump_cm = 40;
cfg_def.minSeqLength = 10;
cfg_def.nMaxNanSkipSequential = 0;
cfg_def.plotOutput = 0;
cfg_def.Qdt = cfg_def.dt/5;
cfg_def.Qboxcar = 5; % if cfg_def.Qdt ~= cfg_def.dt, cfg_def.Qboxcar should be set to cfg_def.dt/cfg_def.Qdt to implement moving window
cfg_def.writeFiles = 1;
cfg_def.removeInterneurons = 0;
cfg_def.keepPosterior = 0;
cfg_def.postCPonly = 0;
cfg_def.use_only_off_track_SWR = 1;

nMaxLaps = 20;
cfg_def.encdecmat = ones(1,nMaxLaps);

cfg = ProcessConfig(cfg_def,cfg_in);

%% load data
please = []; please.load_questionable_cells = cfg.load_questionable_cells;
S = LoadSpikes(please);

LoadExpKeys;
LoadMetadata;

if cfg.removeInterneurons
    cfg_temp = []; cfg_temp.fc = ExpKeys.goodSWR(1);
    csc = LoadCSC(cfg_temp);
    S = RemoveInterneuronsHC([],S,csc);
end

cfg_pos = []; cfg_pos.convFact = ExpKeys.convFact;
pos = LoadPos(cfg_pos); % pos is now in cm
        
%% set up data structs for L, R laps
clear expCond;
expCond(1).label = 'left';
expCond(2).label = 'right';

% match trials
[left,right] = GetMatchedTrials([],metadata,ExpKeys);
expCond(1).t = left;
expCond(2).t = right;

nLapsMax = max(length(expCond(1).t.tstart),length(expCond(2).t.tstart));

% tighter run boundaries: start run based on center photobeam break
evt = getEvents_Tmaze();
for iCond = 1:length(expCond)
    
    this_runStart = expCond(iCond).t.tstart;
    for iT = 1:length(this_runStart)
        next_centerpb_break_idx = nearest_idx3(this_runStart,evt.center_pb,1);
        expCond(iCond).t.tstart = evt.center_pb(next_centerpb_break_idx)'+cfg.trialStartOffset;
    end
end

expCond(1).coord = metadata.coord.coordL_cm;
expCond(2).coord = metadata.coord.coordR_cm;

expCond(1).S = S; % for making tuning curves, "encoding" model
expCond(2).S = S; % this gets restricted by rat running, on track, etc

expCond(1).decS = S; % for decoding
expCond(2).decS = S; % this only gets selections by nSpikes, etc.

expComb.S = S; % combined left and right trials

%% set up output paths
this_fd = pwd;
output_fd = cat(2,pwd,filesep,cfg.output_dir);
base_fn = cat(2,cfg.output_file_prefix,S.cfg.SessionID);

%% linearize paths (snap x,y position samples to nearest point on experimenter-drawn idealized track)
fprintf('Linearizing...');

chp = tsd(0,metadata.coord.chp_cm,{'x','y'}); % make choice point useable by cobebase functions

nCond = length(expCond);
for iCond = 1:nCond
    
%     this_coord.coord = expCond(iCond).coord; this_coord.units = 'px'; this_coord.standardized = 0;
%     expCond(iCond).linpos = LinearizePos([],pos,this_coord);
%     % Compute the coordinate that the choice point corresponds
%     cfg_linpos.outputType = 'dist';
%     
%     chp = tsd(0,metadata.coord.chp_cm,{'x','y'});
%     chp.units = 'cm';
%     chp.cfg.ExpKeys.convFact = ExpKeys.convFact;
%     
%     expCond(iCond).cp = LinearizePos(cfg_linpos,chp,this_coord);
    
    cfg_linpos = []; cfg_linpos.Coord = expCond(iCond).coord;
    expCond(iCond).linpos = LinearizePos_old(cfg_linpos,pos);
    
    % ensure that linpos is now in cm regardless of what we started with
    expCond(iCond).linpos.data = (expCond(iCond).linpos.data ./ length(cfg_linpos.Coord)).*ExpKeys.pathlength;
    
    % get cp in linpos coordinates
    expCond(iCond).cp = LinearizePos_old(cfg_linpos,chp);
    expCond(iCond).cp.data = (expCond(iCond).cp.data ./ length(cfg_linpos.Coord)).*ExpKeys.pathlength;
    
end
    
%% find intervals where rat is running
spd = getLinSpd([],pos); % get speed (in "camera pixels per second")

cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 5;
run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above 5 pix/s

%% exclude positions at beginning and end of track
cfg_track1 = []; cfg_track1.method = 'raw'; cfg_track1.threshold = cfg.trackExcludeStart;
cfg_track2 = []; cfg_track2.method = 'raw'; cfg_track2.operation = '<'; cfg_track2.threshold = ExpKeys.pathlength - cfg.trackExcludeEnd;

both_iv = [];
for iCond = 1:nCond
    track_iv1 = TSDtoIV(cfg_track1,expCond(iCond).linpos);
    if isempty(both_iv)
        both_iv = track_iv1; 
    else
        both_iv = UnionIV([],both_iv,track_iv1);
    end
    
    fh = @(x) restrict(x,track_iv1);
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    track_iv2 = TSDtoIV(cfg_track2,expCond(iCond).linpos);
    
    fh = @(x) restrict(x,track_iv2);
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    both_iv = UnionIV([],both_iv,track_iv2);
end
expComb.S = restrict(expComb.S,both_iv);

%% deal with some rat specific oddities
[~,fn,~] = fileparts(pwd);
switch fn(1:4)
    case 'R042'
        % deal with non-spike sorted chewing intervals
        tf = FindFile('*times.mat');
        if isempty(tf)
            error('No *times.mat file found for R042');
        else
            load(tf);
            chew_iv = iv(t_start*10^-6,t_end*10^-6); % times are in neuralynx timestamps, so convert to s
        end
        
    case 'R044'
        % deal with HS detachments -- remove extended times with no spikes
        cfg_Q = []; cfg_Q.dt = 1;
        Q = MakeQfromS(cfg_Q,S);
        spk_count = tsd(Q.tvec,sum(Q.data));
        
        cfg_det = []; cfg_det.threshold = 0.5; cfg_det.operation = '>'; cfg_det.method = 'raw';
        chew_iv = TSDtoIV(cfg_det,spk_count);
        
    otherwise
        chew_iv = [];
        
end

%% restrict (linearized) position data and spike data to desired intervals
fprintf('Restricting data...');
for iCond = 1:nCond
    
    fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
    expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    
    if ~isempty(chew_iv)
        fh = @(x) restrict(x,chew_iv); % restrict S and linpos to non-detached times
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    end
    
    % also remove cells with insufficient spikes
    [expCond(iCond).S,cell_keep_idx] = removeEmptyCells(expCond(iCond).S);
    expCond(iCond).decS = SelectTS([],expCond(iCond).decS,cell_keep_idx);
    
    spk_count = getSpikeCount([],expCond(iCond).S);
    cell_keep_idx = spk_count >= cfg.minSpikes;
    
    expCond(iCond).S = SelectTS([],expCond(iCond).S,cell_keep_idx);
    expCond(iCond).decS = SelectTS([],expCond(iCond).decS,cell_keep_idx);
end

expComb.S = restrict(expComb.S,run_iv); % don't do trials yet until combined

%% concatenate left and right trials 
temp = expCond(2).linpos;
temp.data = temp.data+ExpKeys.pathlength;
expComb.linpos = UnionTSD([],expCond(1).linpos,temp);
expComb.t = UnionIV([],expCond(1).t,expCond(2).t);
expComb.S = restrict(expComb.S,expComb.t); % restrict with all trials together

% also remove cells with insufficient spikes
[expComb.S,cell_keep_idx] = removeEmptyCells(expComb.S);
expComb.decS = SelectTS([],S,cell_keep_idx);

spk_count = getSpikeCount([],expComb.S);
cell_keep_idx = spk_count >= cfg.minSpikes;

expComb.S = SelectTS([],expComb.S,cell_keep_idx);
expComb.decS = SelectTS([],expComb.decS,cell_keep_idx);

%% get tuning curves & decode
this_TCsmooth = cfg.TCsmooth;

cfg_tc = [];
cfg_tc.binEdges{1} = linspace(0,2*ExpKeys.pathlength,cfg.nBins*2); % ~3cm bins for full maze (smaller for R042...)

if this_TCsmooth ~= 0
    cfg_tc.smoothingKernel = gausskernel(51,this_TCsmooth);
end

enc_S = restrict(expComb.S,expComb.t);
enc_linpos = restrict(expComb.linpos,expComb.t);

expComb.tc = TuningCurves(cfg_tc,enc_S,enc_linpos);
expComb.fields = DetectPlaceCells1D([],expComb.tc.tc); % only used for plotting later

% keep track of choice point
[~,expCond(iCond).cp_bin] = histc(expCond(iCond).cp.data,cfg_tc.binEdges{1});

if cfg.plotOutput
    figure;
    imagesc(expComb.tc.tc);
end

out.n_neurons = size(expComb.tc.tc, 1);
out.tc = expComb.tc;
out.cp_bin = expCond(iCond).cp_bin;

%% Get firing rates for left and right trials before the choice points

trial_nums = length(left.tstart) + length(right.tstart);
avg_FR_per_trial = zeros(length(expComb.S.t), trial_nums);

for iCond = 1:nCond
    cfg_pre_cp = []; cfg_pre_cp.method = 'raw'; cfg_pre_cp.operation = '<';
    cfg_pre_cp.threshold = expCond(iCond).cp.data;

    pre_cp_iv = TSDtoIV(cfg_pre_cp, expCond(iCond).linpos);
    for t_i = 1:length(pre_cp_iv.tstart)
        pre_cp_tstart = pre_cp_iv.tstart(t_i);
        pre_cp_tend = pre_cp_iv.tend(t_i);
        pre_cp_S = restrict(expComb.S, pre_cp_tstart, pre_cp_tend);
        trial_i = (iCond - 1) * (trial_nums / 2) + t_i;
        for n_i = 1:length(expComb.S.t)
            avg_FR_per_trial(n_i, trial_i) = length(pre_cp_S.t{n_i}) / (pre_cp_tend - pre_cp_tstart);
        end
    end
end

% figure; imagesc(avg_FR_per_trial); colorbar;
% xlabel('trial'); ylabel('neurons')

L_FR = mean(avg_FR_per_trial(:, 1:trial_nums/2), 2);
R_FR = mean(avg_FR_per_trial(:, trial_nums/2+1:end), 2);
FR_diff = abs(L_FR - R_FR);
out.FR_diff = FR_diff;
% figure; imagesc([L_FR, R_FR, FR_diff]); colorbar;

%% Comparing FR of left and right with shuffles

FR_n_shuffles = 1000;
shuf_FR_diff = zeros(length(expComb.S.t), FR_n_shuffles);

for shuf_i = 1:FR_n_shuffles
    rand_indices = randperm(trial_nums);
    shuf_FR_per_trial = avg_FR_per_trial(:, rand_indices);
    s_L_FR = mean(shuf_FR_per_trial(:, 1:trial_nums/2), 2);
    s_R_FR = mean(shuf_FR_per_trial(:, trial_nums/2+1:end), 2);
    shuf_FR_diff(:, shuf_i) = abs(s_L_FR - s_R_FR);
end

% figure; imagesc([shuf_FR_diff]); colorbar;

FR_diff_zscore = (FR_diff - mean(shuf_FR_diff, 2, 'omitnan')) ./ std(shuf_FR_diff, 0, 2, 'omitnan');
out.FR_diff_zscore = FR_diff_zscore;
% figure; scatter(FR_diff, FR_diff_zscore, 'filled');
% xlabel('FR diff'); ylabel('FR diff zscore');
% title(['FR diff zscore: ', num2str(median(abs(out.FR_diff_zscore), 'omitnan'))]);

%% Q-mat
cfg_Q = [];
cfg_Q.dt = cfg.Qdt;
cfg_Q.boxcar_size = cfg.Qboxcar;

if cfg.QsmoothSD == 0
    cfg_Q.smooth = [];
else
    cfg_Q.smooth = 'gauss';
    cfg_Q.gausswin_sd = cfg.QsmoothSD;
end

expComb.Q = MakeQfromS(cfg_Q,expComb.decS);

% raw Q-mat, to get nNeurons later
cfg_Qraw = []; cfg_Qraw.smooth = []; cfg_Qraw.dt = cfg_Q.dt; cfg_Qraw.boxcar_size = cfg_Q.boxcar_size;
expComb.Qraw = MakeQfromS(cfg_Qraw,expComb.decS);
expComb.nNeurons = sum(expComb.Qraw.data >= 1);

%% decode
clear csc pos S;
cfg_decode = [];
cfg_decode.nMinSpikes = cfg.dt;
cfg_decode.excludeMethod = 'frate';
expComb.P = DecodeZ(cfg_decode,expComb.Q,expComb.tc.tc); % full decoded probability distribution

%% get MAP & plot output
[~,map] = max(expComb.P.data);
toss_idx = isnan(nansum(expComb.P.data));
map(toss_idx) = NaN;

if cfg.plotOutput
    imagesc(expComb.P.tvec,1:size(expComb.P.data,1),expComb.P.data);
    hold on;
    plot(expComb.P.tvec,map,'.w');
    plot(expComb.linpos.tvec,expComb.tc.usr.pos_idx,'og');
    
    figure;
    subplot(221);
    hist(map,cfg.nBins*2);
    title('MAP histogram');
    
    subplot(222);
    hist(nansum(expComb.P.data));
    title('summed posterior histogram');

end


%% quantify decoding accuracy on RUN
this_trueZ = tsd(expComb.linpos.tvec,expComb.tc.usr.pos_idx); % true position in units of bins (as established by tuning curves)
cfg_err = []; cfg_err.mode = 'max';

keep_idx = unique(nearest_idx3(this_trueZ.tvec,expComb.P.tvec)); % match up decoding with true positions
this_Pscore = expComb.P;
this_Pscore.tvec = this_Pscore.tvec(keep_idx);
this_Pscore.data = this_Pscore.data(:,keep_idx);
[expComb.Perr,expComb.confMat] = DecodeErrorZ(cfg_err,this_Pscore,this_trueZ);

if cfg.plotOutput
    imagesc(expComb.confMat.full);
    caxis([0 0.1]);

    hold on;
    h = plot([112 112],[1 223],'w');
    h2 = plot([1 223],[112 112],'w');
end

%% compare quadrants
cfy = @(x) x(:);
this_mat = expComb.confMat.thr;

keep = ~isnan(nansum(expComb.confMat.thr')); % find non-nan columns
left_pts = [1 expCond(iCond).cp_bin cfg.nBins-1]; np_l1 = sum(keep(left_pts(1):left_pts(2))); np_l2 = sum(keep(left_pts(2)+1:left_pts(3)));
right_pts = [cfg.nBins+left_pts]; np_r1 = sum(keep(right_pts(1):right_pts(2))); np_r2 = sum(keep(right_pts(2)+1:right_pts(3)));

out.left_preCP_correct = nansum(cfy(this_mat(left_pts(1):left_pts(2),left_pts(1):left_pts(2))))./np_l1;
out.left_preCP_incorrect = nansum(cfy(this_mat(left_pts(1):left_pts(2),right_pts(1):right_pts(2))))./np_l1; % points classified as OPPOSITE trial equivalent
out.right_preCP_correct = nansum(cfy(this_mat(right_pts(1):right_pts(2),right_pts(1):right_pts(2))))./np_r1;
out.right_preCP_incorrect = nansum(cfy(this_mat(right_pts(1):right_pts(2),left_pts(1):left_pts(2))))./np_r1;
fprintf('PreCP: left+ %.3f, left- %.3f, right+ %.3f, right- %.3f\n',out.left_preCP_correct,out.left_preCP_incorrect,out.right_preCP_correct,out.right_preCP_incorrect);

out.left_postCP_correct = nansum(cfy(this_mat(left_pts(2)+1:left_pts(3),left_pts(2)+1:left_pts(3))))./np_l2;
out.left_postCP_incorrect = nansum(cfy(this_mat(left_pts(2)+1:left_pts(3),right_pts(2)+1:right_pts(3))))./np_l2;
out.right_postCP_correct = nansum(cfy(this_mat(right_pts(2)+1:right_pts(3),right_pts(2)+1:right_pts(3))))./np_r2;
out.right_postCP_incorrect = nansum(cfy(this_mat(right_pts(2)+1:right_pts(3),left_pts(2)+1:left_pts(3))))./np_r2;
fprintf('PostCP: left+ %.3f, left- %.3f, right+ %.3f, right- %.3f\n',out.left_postCP_correct,out.left_postCP_incorrect,out.right_postCP_correct,out.right_postCP_incorrect);

out.left_chance_pre = np_l1./(np_l1+np_l2+np_r1+np_r2); out.left_chance_post = np_l2./(np_l1+np_l2+np_r1+np_r2); % analytical chance levels based on number of bins
out.right_chance_pre = np_r1./(np_l1+np_l2+np_r1+np_r2); out.right_chance_post = np_r2./(np_l1+np_l2+np_r1+np_r2); % analytical chance levels based on number of bins

%% now decode SWR vectors
%
cfg_shuf_LR.dt = 0.05; % this sets the binsize in decoding
LoadCandidates;

if cfg.use_only_off_track_SWR
    rest_iv = metadata.taskvars.rest_iv;
    off_track_iv = rest_iv;
    off_track_iv.tstart = [ExpKeys.prerecord(1); rest_iv.tstart; ExpKeys.postrecord(1)];
    off_track_iv.tend = [ExpKeys.prerecord(2); rest_iv.tend; ExpKeys.postrecord(2)];
    evt = restrict(evt, off_track_iv);
end

% make Q-matrix based on events
cfg_Q = [];
%cfg_Q.tvec_edges = sort(cat(1,evt.tstart,evt.tend)); % raw events, have variable length...
tvec_centers = IVcenters(evt); % alternative: symmetric around center
cfg_Q.tvec_edges = sort(cat(1,tvec_centers-cfg_shuf_LR.dt/2,tvec_centers+cfg_shuf_LR.dt/2));

Q_SWR = MakeQfromS(cfg_Q,expComb.decS);
Q_SWR.data = Q_SWR.data(:,1:2:end);
Q_SWR.tvec = Q_SWR.tvec(1:2:end);

% need to do some selection on minimum number of cells?
nActiveCells = sum(Q_SWR.data > 0);
keep = nActiveCells >= cfg.nMinNeurons;
Q_SWR.tvec = Q_SWR.tvec(keep); Q_SWR.data = Q_SWR.data(:,keep);
% keep raw tvecs for subsequent analysis
Q_SWR.raw_tvec = Q_SWR.tvec;

Q_SWR.tvec = cfg_shuf_LR.dt*(1:length(Q_SWR.tvec)); 
tvec_centers = tvec_centers(keep);

% decode
cfg_decode = [];
cfg_decode.nMinSpikes = cfg.dt;
cfg_decode.excludeMethod = 'frate';
expComb.P_SWR = DecodeZ(cfg_decode,Q_SWR,expComb.tc.tc); % full decoded probability distribution

% obtain log odds
div = ceil(size(expComb.tc.tc,2)/2); % divider between L & R tuning curves
if cfg.postCPonly
    this_cp = expCond(2).cp_bin;
    pL = nansum(expComb.P_SWR.data(this_cp+1:div-1,:)); pR = nansum(expComb.P_SWR.data(div+this_cp+1:end,:));
else
    pL = nansum(expComb.P_SWR.data(1:div-1,:)); pR = nansum(expComb.P_SWR.data(div:end,:));
end
actual_odds = log2(pL./pR); % this measure is not ideal because of Infs, but pL-pR gives similar results
actual_diff = pL - pR;

%%
out.actual_pL = pL; out.actual_pR = pR;
out.raw_tvec = Q_SWR.raw_tvec;

%% L vs R shuffle
cfg_shuf_LR = [];
cfg_shuf_LR.nShuffles = 1000;

nCells = size(expComb.tc.tc,1);

clear shuf_odds;
clear shuf_diff;
for iShuf = cfg_shuf_LR.nShuffles:-1:1
   
    this_tc = expComb.tc.tc;
    
    % figure out which TCs to swap
    swap = logical(randi(2,nCells,1)-1);
    this_tc(swap,1:div) = expComb.tc.tc(swap,div:end);
    this_tc(swap,div:end) = expComb.tc.tc(swap,1:div);
    
    expComb.P_SWR = DecodeZ(cfg_decode,Q_SWR,this_tc); % full decoded probability distribution
    
    % obtain log odds
    div = ceil(size(expComb.tc.tc,2)/2); % divider between L & R tuning curves
    if cfg.postCPonly
        pL = nansum(expComb.P_SWR.data(this_cp+1:div-1,:)); pR = nansum(expComb.P_SWR.data(div+this_cp+1:end,:));
    else
        pL = nansum(expComb.P_SWR.data(1:div-1,:)); pR = nansum(expComb.P_SWR.data(div:end,:));
    end
    shuf_odds(iShuf,:) = log2(pL./pR);
    shuf_diff(iShuf,:) = pL - pR;
end
out.shuf_perc_odds = nansum(shuf_odds < repmat(actual_odds,[cfg_shuf_LR.nShuffles 1]))./cfg_shuf_LR.nShuffles;
out.shuf_perc_odds(isnan(actual_odds)) = NaN;
out.shuf_z_odds = (actual_odds-nanmean(shuf_odds))./nanstd(shuf_odds);
out.shuf_z_odds(isnan(actual_odds)) = NaN;

out.shuf_perc_diff = nansum(shuf_diff < repmat(actual_diff,[cfg_shuf_LR.nShuffles 1]))./cfg_shuf_LR.nShuffles;
out.shuf_perc_diff(isnan(actual_diff)) = NaN;
out.shuf_z_diff = (actual_diff-nanmean(shuf_diff))./nanstd(shuf_diff);
out.shuf_z_diff(isnan(actual_diff)) = NaN;

out.tvec = tvec_centers;
out.behav_sequence = metadata.taskvars.sequence;
% odds are left over right. So percentile = big means left, percentile =
% small means right.

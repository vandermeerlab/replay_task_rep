function out = Get_Gupta_DecSeqCombined(cfg_in, TC)

%% General configs
cfg_def = [];
cfg_def.load_questionable_cells = 1;
cfg_def.removeInterneurons = 1;
cfg_def.minSpikes = 25;
cfg_def.nMinNeurons = 0;
cfg_def.dt = 0.1;
cfg_def.Qdt = cfg_def.dt/5;
cfg_def.Qboxcar = 5;
cfg_def.postCPonly = 0;

cfg = ProcessConfig(cfg_def,cfg_in);
%% load spike data
LoadMetadata;
LoadExpKeys();

please = []; please.load_questionable_cells = cfg.load_questionable_cells;
S = LoadSpikes(please);

if cfg.removeInterneurons
    cfg_lfp = []; cfg_lfp.fc = ExpKeys.goodSWR(1);
    lfp = LoadCSC(cfg_lfp);
    S = RemoveInterneuronsHC([], S, lfp);
end

decS = S; % for decoding, this only gets selections by nSpikes, etc.

% remove cells with insufficient spikes
[~, cell_keep_idx] = removeEmptyCells(S);
decS = SelectTS([], decS, cell_keep_idx);

spk_count = getSpikeCount([], S);
cell_keep_idx = spk_count >= cfg.minSpikes;

decS = SelectTS([], decS, cell_keep_idx);

%% record relevant task info
out.contigency = ExpKeys.maze;
out.switch_time = metadata.SwitchTime;

switch_indices = find(metadata.taskvars.trial_iv.tstart > out.switch_time);
switch_idx = switch_indices(1);
out.switch_idx = switch_idx;

out.behav_sequence = metadata.taskvars.sequence;
out.trial_iv = metadata.taskvars.trial_iv;
out.TimeOffTrack = metadata.TimeOffTrack;

%% Get firing rates for left and right trials before the choice points
L_nums = length(metadata.taskvars.trial_iv_L.tstart);
R_nums = length(metadata.taskvars.trial_iv_R.tstart);
trial_nums = L_nums + R_nums;
avg_FR_per_trial = zeros(length(TC.combined.S.t), trial_nums);

all_TCs = {TC.left, TC.right};
for iCond = 1:length(all_TCs)
    cfg_pre_cp = []; cfg_pre_cp.method = 'raw'; cfg_pre_cp.operation = '<';
    cfg_pre_cp.threshold = all_TCs{iCond}.cp.data;

    pre_cp_iv = TSDtoIV(cfg_pre_cp, all_TCs{iCond}.linpos);
    for t_i = 1:length(pre_cp_iv.tstart)
        pre_cp_tstart = pre_cp_iv.tstart(t_i);
        pre_cp_tend = pre_cp_iv.tend(t_i);
        pre_cp_S = restrict(TC.combined.S, pre_cp_tstart, pre_cp_tend);
        trial_i = (iCond - 1) * L_nums + t_i;
        for n_i = 1:length(TC.combined.S.t)
            avg_FR_per_trial(n_i, trial_i) = length(pre_cp_S.t{n_i}) / (pre_cp_tend - pre_cp_tstart);
        end
    end
end

% figure; imagesc(avg_FR_per_trial); colorbar;
% xlabel('trial'); ylabel('neurons')

L_FR = mean(avg_FR_per_trial(:, 1:L_nums), 2);
R_FR = mean(avg_FR_per_trial(:, L_nums+1:end), 2);
FR_diff = abs(L_FR - R_FR);
out.FR_diff = FR_diff;
% figure; imagesc([L_FR, R_FR, FR_diff]); colorbar;

%% Q-mat
cfg_Q = [];
cfg_Q.dt = cfg.Qdt;
cfg_Q.boxcar_size = cfg.Qboxcar;
cfg_Q.smooth = [];

Q = MakeQfromS(cfg_Q, decS);

%% decode
clear csc pos S;
cfg_decode = [];
cfg_decode.nMinSpikes = cfg.dt;
cfg_decode.excludeMethod = 'frate';
P = DecodeZ(cfg_decode,Q,TC.combined.tc.tc); % full decoded probability distribution

%% quantify decoding accuracy on RUN
this_trueZ = tsd(TC.combined.linpos.tvec,TC.combined.tc.usr.pos_idx); % true position in units of bins (as established by tuning curves)
cfg_err = []; cfg_err.mode = 'max';

keep_idx = unique(nearest_idx3(this_trueZ.tvec,P.tvec)); % match up decoding with true positions
this_Pscore = P;
this_Pscore.tvec = this_Pscore.tvec(keep_idx);
this_Pscore.data = this_Pscore.data(:,keep_idx);
[Perr, confMat] = DecodeErrorZ(cfg_err,this_Pscore,this_trueZ);

%% compare quadrants
cfy = @(x) x(:);
this_mat = confMat.thr;

keep = ~isnan(nansum(confMat.thr')); % find non-nan columns
left_pts = [1 TC.left.cp_bin TC.combined.nBins/2]; np_l1 = sum(keep(left_pts(1):left_pts(2))); np_l2 = sum(keep(left_pts(2)+1:left_pts(3)));
right_pts = [TC.combined.nBins/2+1 TC.combined.nBins/2+TC.right.cp_bin, TC.combined.nBins]; np_r1 = sum(keep(right_pts(1):right_pts(2))); np_r2 = sum(keep(right_pts(2)+1:right_pts(3)));

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

% make Q-matrix based on events
cfg_Q = [];
%cfg_Q.tvec_edges = sort(cat(1,evt.tstart,evt.tend)); % raw events, have variable length...
tvec_centers = IVcenters(evt); % alternative: symmetric around center
cfg_Q.tvec_edges = sort(cat(1,tvec_centers-cfg_shuf_LR.dt/2,tvec_centers+cfg_shuf_LR.dt/2));

Q_SWR = MakeQfromS(cfg_Q,decS);
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
P_SWR = DecodeZ(cfg_decode,Q_SWR,TC.combined.tc.tc); % full decoded probability distribution

% obtain log odds
div = ceil(size(TC.combined.tc.tc,2)/2); % divider between L & R tuning curves
if cfg.postCPonly
    this_cp = min(TC.left.cp_bin, TC.right.cp_bin);
    pL = nansum(P_SWR.data(this_cp+1:div-1,:)); pR = nansum(P_SWR.data(div+this_cp+1:end,:));
else
    pL = nansum(P_SWR.data(1:div-1,:)); pR = nansum(P_SWR.data(div:end,:));
end
actual_odds = log2(pL./pR); % this measure is not ideal because of Infs, but pL-pR gives similar results
actual_diff = pL - pR;

%%
out.actual_pL = pL; out.actual_pR = pR;
out.actual_odds = actual_odds; out.actual_diff = actual_diff;
out.raw_tvec = Q_SWR.raw_tvec;

%% L vs R shuffle
cfg_shuf_LR = [];
cfg_shuf_LR.nShuffles = 1000;

nCells = size(TC.combined.tc.tc,1);

clear shuf_odds;
clear shuf_diff;
for iShuf = cfg_shuf_LR.nShuffles:-1:1

    this_tc = TC.combined.tc.tc;

    % figure out which TCs to swap
    swap = logical(randi(2,nCells,1)-1);
    this_tc(swap,1:div) = TC.combined.tc.tc(swap,div+1:end);
    this_tc(swap,div+1:end) = TC.combined.tc.tc(swap,1:div);

    sf_P_SWR = DecodeZ(cfg_decode,Q_SWR,this_tc); % full decoded probability distribution

    % obtain log odds
    div = ceil(size(TC.combined.tc.tc,2)/2); % divider between L & R tuning curves
    if cfg.postCPonly
        pL = nansum(sf_P_SWR.data(this_cp+1:div-1,:)); pR = nansum(sf_P_SWR.data(div+this_cp+1:end,:));
    else
        pL = nansum(sf_P_SWR.data(1:div-1,:)); pR = nansum(sf_P_SWR.data(div:end,:));
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
% odds are left over right. So percentile = big means left, percentile =
% small means right.

function out = Get_Gupta_DecSeqCombined(cfg_in, TC)

%% General configs
cfg_def = [];
cfg_def.load_questionable_cells = 1;
cfg_def.removeInterneurons = 1;
cfg_def.minSpikes = 25;
cfg_def.nMinNeurons = 0;
cfg_def.dt = 0.1;
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
    this_cp = expCond(2).cp_bin;
    pL = nansum(P_SWR.data(this_cp+1:div-1,:)); pR = nansum(P_SWR.data(div+this_cp+1:end,:));
else
    pL = nansum(P_SWR.data(1:div-1,:)); pR = nansum(P_SWR.data(div:end,:));
end
actual_odds = log2(pL./pR); % this measure is not ideal because of Infs, but pL-pR gives similar results
actual_diff = pL - pR;

%%
out.actual_pL = pL; out.actual_pR = pR;
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

    expComb.P_SWR = DecodeZ(cfg_decode,Q_SWR,this_tc); % full decoded probability distribution

    % obtain log odds
    div = ceil(size(TC.combined.tc.tc,2)/2); % divider between L & R tuning curves
    if cfg.postCPonly
        pL = nansum(P_SWR.data(this_cp+1:div-1,:)); pR = nansum(P_SWR.data(div+this_cp+1:end,:));
    else
        pL = nansum(P_SWR.data(1:div-1,:)); pR = nansum(P_SWR.data(div:end,:));
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

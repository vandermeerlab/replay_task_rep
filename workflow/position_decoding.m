%% General configs
cfg.load_questionable_cells = 1;
cfg.minSpikes = 25;
cfg.dt = 0.1;
cfg.Qdt = cfg.dt/5;
cfg.Qboxcar = 5; % if cfg_def.Qdt ~= cfg_def.dt, cfg_def.Qboxcar should be set to cfg_def.dt/cfg_def.Qdt to implement moving window
cfg.QsmoothSD = 0;
cfg.plotOutput = 1;

%% load spike data
please = []; please.load_questionable_cells = cfg.load_questionable_cells;
S = LoadSpikes(please);

decS = S; % for decoding, this only gets selections by nSpikes, etc.

% remove cells with insufficient spikes
[S, cell_keep_idx] = removeEmptyCells(S);
decS = SelectTS([], decS, cell_keep_idx);

spk_count = getSpikeCount([], S);
cell_keep_idx = spk_count >= cfg.minSpikes;

S = SelectTS([], S, cell_keep_idx);
decS = SelectTS([], decS, cell_keep_idx);

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

Q = MakeQfromS(cfg_Q, decS);

%% decode
clear S;
cfg_decode = [];
cfg_decode.nMinSpikes = cfg.dt;
cfg_decode.excludeMethod = 'frate';
P = DecodeZ(cfg_decode, Q, TC.combined.tc.tc); % full decoded probability distribution

%% get MAP & plot output
[~,map] = max(P.data);
toss_idx = isnan(nansum(P.data));
map(toss_idx) = NaN;

if cfg.plotOutput
    imagesc(P.tvec,1:size(P.data,1),P.data);
    hold on;
    plot(P.tvec,map,'.w');
    plot(TC.combined.linpos.tvec, TC.combined.tc.usr.pos_idx,'og');
    
    figure;
    subplot(221);
    hist(map,TC.combined.nBins);
    title('MAP histogram');
    
    subplot(222);
    hist(nansum(P.data));
    title('summed posterior histogram');

end

%% quantify decoding accuracy on RUN
trueZ = tsd(TC.combined.linpos.tvec, TC.combined.tc.usr.pos_idx); % true position in units of bins (as established by tuning curves)
cfg_err = []; cfg_err.mode = 'max';

keep_idx = unique(nearest_idx3(trueZ.tvec, P.tvec)); % match up decoding with true positions
Pscore = P;
Pscore.tvec = Pscore.tvec(keep_idx);
Pscore.data = Pscore.data(:,keep_idx);
[Perr, confMat] = DecodeErrorZ(cfg_err, Pscore, trueZ);

if cfg.plotOutput
    figure;
    imagesc(confMat.full);
    caxis([0 0.1]);
    
    hold on;
    h = plot([100 100],[1 200],'w');
    h2 = plot([1 200],[100 100],'w');
    colorbar;
    ylabel('actual bins')
    xlabel('decoded bins')
    set(gca,'FontSize',12)
end
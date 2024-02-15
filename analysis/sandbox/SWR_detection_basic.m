%% some hippocampus data
cfg = []; cfg.fc = ExpKeys.goodSWR;
lfp = LoadCSC(cfg);
 
%% filter in SWR band
cfg = [];
cfg.f = [140 220];
cfg.display_filter = 0;
 
SWRf = FilterLFP(cfg,lfp);
 
%% obtain power and z-score it
SWRp = LFPpower([],SWRf);
SWRp_z = zscore_tsd(SWRp);
 
%% detect events
cfg = [];
cfg.method = 'raw';
cfg.threshold = 3;
cfg.operation =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.minlen = 0.05; % minimum interval length
 
SWR_evt = TSDtoIV(cfg,SWRp_z);
 
%% to each event, add a field with the max z-scored power (for later selection)
cfg = [];
cfg.method = 'max'; % 'min', 'mean'
cfg.label = 'maxSWRp'; % what to call this in iv, i.e. usr.label
 
SWR_evt = AddTSDtoIV(cfg,SWR_evt,SWRp_z);
 
%% select only those events of >5 z-scored power
cfg = [];
cfg.operation = '>';
cfg.threshold = 5;
 
SWR_evt = SelectIV(cfg,SWR_evt,'maxSWRp');
 
%% plot events in highlighted on top of full lfp
PlotTSDfromIV([],SWR_evt,lfp);

%%
cfg = [];
cfg.display = 'iv';
cfg.mode = 'center';
cfg.fgcol = 'k';
 
PlotTSDfromIV(cfg,SWR_evt,lfp);
%% ..hold on (highlight edges of event on top of previous plot)
cfg = [];
cfg.display = 'iv';
cfg.fgcol = 'r';
 
PlotTSDfromIV(cfg,SWR_evt,lfp);

%% restrict the data to interval of interest
this_iv = iv(SWR_evt.tstart, SWR_evt.tend);
S_r = restrict(S,this_iv);
csc_r = restrict(csc,this_iv);
 
%% plot spikes
SET_spkY = 0.4; % parameter to set spike height
 
figure; hold on;
for iC = 1:length(S_r.t)
    if ~isempty(S_r.t{iC})
        plot([S_r.t{iC} S_r.t{iC}],[iC-SET_spkY iC+SET_spkY],'k');
    end
end % of cells
ylabel('neuron #');
 
%% add a LFP
SET_cscY = [-5 0];
plot(csc_r.tvec,rescale(csc_r.data,SET_cscY(1),SET_cscY(2)),'r');
set(gca,'YTickLabel',{});
 
%% add multi-unit activity in separate axes
ax1 = gca; ax1_pos = get(ax1,'Position');
ax2 = axes('Position',ax1_pos); % new axes with same position as first
 
cfg = []; cfg.tvec = csc.tvec; cfg.sigma = 0.1;
mua = getMUA(cfg,S); % obtain multi-unit activity
 
xr = get(ax1,'XLim');
mua_r = restrict(mua,xr(1),xr(2)); % only keep the data we need
 
axes(ax2); % set current axes to the second set, so what follows happens there
mua_hdl = plot(mua_r.tvec,mua_r.data,'Color',[0.7 0.7 0.7]);
 
set(gca,'YAxisLocation','right','Box','off','XTick',[],'Color','none','YColor',[0.5 0.5 0.5])
ylabel('multi-unit activity (spk/s)');
linkaxes([ax1 ax2],'x'); 
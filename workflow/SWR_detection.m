%%
% Do you want to load the 5's?
cfg.load_questionable_cells = 1;

% Which SWR detector do you want to use? 
cfg.SWRmethod = 'HT'; % 'AM' for amSWR (frequency content similarity), 'HT' for OldWizard (hilbert transform), 'TR' for photonic (transient detection), or 'none' to skip

% Which MUA detector do you want to use?
cfg.MUAmethod = 'none'; % 'AM' for amMUA, or 'none' for skip MUA detection

cfg.weightby = 'amplitude'; % this applies to 'AM' amSWR and 'TR' photonic, but not 'HT' OldWizard

% If using amSWR, how fast do you want it to go?
cfg.stepSize = 4; % this applies to 'AM' SWR method only

% Threshold for making intervals
cfg.DetectorThreshold = 3; % the threshold you want for generating IV data

% How is the threshold done? 'raw','zscore'
cfg.ThreshMethod = 'zscore';

% Minimum duration of intervals
cfg.mindur = 0.02; % in seconds

% Disclude events when the rat was moving faster than this
cfg.SpeedLimit = 10; % ( suggested 10 ) pixels per second, if [] no speed thresholding

% Disclude events with high theta power
cfg.ThetaThreshold = 2; % ( suggested 2 ) power std above mean, if [] no theta limiting

% Minimum number of active cells for the event to be kept
cfg.minCells = 5;

% Amount to add to the interval (catch borderline missed spikes)
cfg.expandIV = [0 0]; % in seconds (see ResizeIV)

cfg.allowOverlap = 0; % don't allow the expanded intervals to overlap one another

%%
SWR_evt = GenCandidateEvents(cfg);
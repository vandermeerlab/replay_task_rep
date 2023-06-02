%%
data_paths = getGuptaDataPath([]);

%%
cfg_tc = [];
cfg_tc.use_Gupta_data = 1;
cfg_tc.use_matched_trials = 0;
cfg_tc.removeInterneurons = 1;
TC = get_tuning_curve(cfg_tc, data_paths{1});
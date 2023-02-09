%%
data_paths = getGuptaDataPath([]);
for p_i = 1:length(data_paths)
    cd(data_paths{p_i});
end

%%
cfg_tc = [];
cfg_tc.use_Gupta_data = 1;
cfg_def.use_matched_trials = 0;
TC = get_tuning_curve(cfg_tc, data_paths{1});
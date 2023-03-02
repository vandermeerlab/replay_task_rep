%%
data_paths = getGuptaDataPath([]);
for p_i = 1:length(data_paths)
    cd(data_paths{p_i});
end

%%
cfg_tc = [];
cfg_tc.use_Gupta_data = 1;
cfg_tc.use_matched_trials = 0;
cfg_tc.removeInterneurons = 0;
enc_TC = get_tuning_curve(cfg_tc, data_paths{1});
TC.left.tc = enc_TC.left.tc.tc;
TC.right.tc = enc_TC.right.tc.tc;
TC.combined.tc = [TC.left.tc, TC.right.tc];
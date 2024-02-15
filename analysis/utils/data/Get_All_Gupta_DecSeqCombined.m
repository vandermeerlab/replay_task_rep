% batch script to generate non-sequence decoding output for each session
%
% the only configuration setting you may want to change here is:
%
% cfg_tc.removeInterneurons: set to 0 (default) to include all cells for
% decoding, set to 1 to remove putative interneurons (> 5 Hz firing rate)

%%
data_paths = getGuptaDataPath([]);

%%
for p_i = 1:length(data_paths)
    cd(data_paths{p_i});
    cfg_tc = [];
    cfg_tc.use_Gupta_data = 1;
    cfg_tc.use_matched_trials = 0;
    cfg_tc.removeInterneurons = 1;
    TC = get_tuning_curve(cfg_tc, data_paths{p_i});

    cfg_decSeq = [];
    cfg_decSeq.postCPonly = 0; % exclude central stem of T-maze up until the choice point?
    out{p_i} = Get_Gupta_DecSeqCombined(cfg_decSeq, TC);
end

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('~~~                      End of script run                          ~~~')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
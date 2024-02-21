%% Calculate the distribution of time difference as a function of delta events
devt_max = 10;
sig_events_only = 1;

SWR_t_diff_by_devt = cell(1, length(out));

for sess_i = 1:length(out)
    SWR_data = out{sess_i}.actual_pL - out{sess_i}.actual_pR;
    SWR_times = out{sess_i}.tvec;

    if sig_events_only
        L_sig_odds = out{sess_i}.shuf_perc_odds >= 0.95;
        R_sig_odds = out{sess_i}.shuf_perc_odds <= 0.05;
        SWR_data = SWR_data(L_sig_odds | R_sig_odds);
        SWR_times = SWR_times(L_sig_odds | R_sig_odds);
    end
    SWR_times = SWR_times(~isnan(SWR_data));

    cfg_time = [];
    cfg_time.devt_max = devt_max;
    [SWR_t_diff_by_devt{sess_i}] = get_time_dist_delta_events(cfg_time, SWR_times);
end

%% Plotting the distribution for each session
for sess_i = 1:length(out)
    SWR_t_diffs = [];
    devt_labels = {};
    devt_order = cell(devt_max, 1);
    
    for d_i = 1:devt_max
        SWR_t_diffs = [SWR_t_diffs; SWR_t_diff_by_devt{sess_i}{d_i}];
        for t_i = 1:length(SWR_t_diff_by_devt{sess_i}{d_i})
            devt_labels{end+1} = num2str(d_i);
        end
        devt_order{d_i} = num2str(d_i);
    end
    % Violin plot
    figure
    vs = violinplot(log10(SWR_t_diffs), devt_labels', 'GroupOrder', devt_order);
    xlabel('delta events');
    ylabel('log10 time elapsed between events(s)');
    ylim([-1.5, 3]);
end

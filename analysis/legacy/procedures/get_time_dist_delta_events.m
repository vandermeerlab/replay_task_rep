function [SWR_t_diff_by_devt] = get_time_dist_delta_events(cfg_in, SWR_times)

    % Calculate the distribution of time difference as a function of delta events
    cfg_def = [];
    cfg_def.devt_max = 10;
    cfg = ProcessConfig(cfg_def,cfg_in);
    
    SWR_t_diff_by_devt = cell(cfg.devt_max, 1);
    for d_i = 1:cfg.devt_max
        for evt_i = 1:length(SWR_times)-1
            if evt_i + d_i <= length(SWR_times)
                SWR_t_diff = SWR_times(evt_i + d_i) - SWR_times(evt_i);
                if SWR_t_diff > 0.1
                    SWR_t_diff_by_devt{d_i} = [SWR_t_diff_by_devt{d_i}; SWR_t_diff];
                end
            end
        end
    end

end
function [p_switch, num_events_bin] = calculate_p_switch_by_time(cfg_in, SWR_data, SWR_times)

cfg_def = [];
cfg_def.bin_egdes = -1.2:0.2:1.8;

cfg = ProcessConfig(cfg_def,cfg_in);

SWR_t_diffs = [];
SWR_t_switch = [];

for i = 1:length(SWR_times)-1
    %     for j = i+1:length(sig_SWR_idx)
    %         if i ~= j
    SWR_t_diff = SWR_times(i+1) - SWR_times(i);
    if SWR_t_diff > 0.1
        SWR_t_diffs = [SWR_t_diffs, SWR_t_diff];
        L_R_dt2 = SWR_data(i+1);
        L_R_dt1 = SWR_data(i);
        % Determine stick or switch
        if L_R_dt1 > 0 && L_R_dt2 > 0
            % Both are L
            SWR_t_switch = [SWR_t_switch, 0];
        elseif L_R_dt1 < 0 && L_R_dt2 < 0
            % Both are R
            SWR_t_switch = [SWR_t_switch, 0];
        else
            SWR_t_switch = [SWR_t_switch, 1];
        end
    end
    %         end
    %     end
end

p_switch = zeros(1, length(cfg.bin_egdes)-1);
num_events_bin = zeros(1, length(cfg.bin_egdes)-1);

t_diffs_bin_indices = discretize(log10(SWR_t_diffs), cfg.bin_egdes);
for bin_i = 1:length(p_switch)
    p_switch_bin = SWR_t_switch(t_diffs_bin_indices == bin_i);
    if isempty(p_switch_bin)
        p_switch(bin_i) = NaN;
        num_events_bin(bin_i) = 0;
    else
        p_switch(bin_i) = sum(p_switch_bin) / length(p_switch_bin);
        num_events_bin(bin_i) = length(p_switch_bin);
    end
end

end
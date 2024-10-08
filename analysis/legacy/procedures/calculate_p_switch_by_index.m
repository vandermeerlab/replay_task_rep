function [p_switch] = calculate_p_switch_by_index(cfg_in, SWR_data)

% Calculate the probablity of switching across delta event index
% Takes SWR_data, size of [1, n_SWR_events] whose element is L if > 0, R if
% < 0
% Return p_switch, size of [1, devt_max] whose elements containing the
% proportion of switching between events across delta indices

cfg_def = [];
cfg_def.devt_max = 50;
% cfg_def.devt_nbin = 1;

cfg = ProcessConfig(cfg_def,cfg_in);

for d_i = 1:cfg.devt_max
    p_switch_devt = [];
    for evt_i = 1:length(SWR_data)-1
        if evt_i + d_i <= length(SWR_data)
            SWR_evt1 = SWR_data(evt_i);
            SWR_evt2 = SWR_data(evt_i + d_i);
            if ~isnan(SWR_evt1) && ~isnan(SWR_evt2)
                % Both events are L
                if SWR_evt1 > 0 && SWR_evt2 > 0
                    p_switch_devt = [p_switch_devt, 0];
                % Both events are R
                elseif SWR_evt1 < 0 && SWR_evt2 < 0
                    p_switch_devt = [p_switch_devt, 0];
                % Switch (R then L)
                elseif SWR_evt1 < 0 && SWR_evt2 > 0
                    p_switch_devt = [p_switch_devt, 1];
                % Switch (L then R)
                elseif SWR_evt1 > 0 && SWR_evt2 < 0
                    p_switch_devt = [p_switch_devt, 1];
                end
            end
        end
        p_switch(d_i) = sum(p_switch_devt) / length(p_switch_devt);
        n_delta_events(d_i) = length(p_switch_devt);
    end
end

% if cfg.devt_nbin > 1
%     n_bins = cfg.devt_max / cfg.devt_nbin;
%     p_switch_bin = zeros(1, n_bins);
%     for bin_i = 0:n_bins-1
%         devt_i_start = bin_i * cfg.devt_nbin + 1;
%         devt_i_end = devt_i_start + cfg.devt_nbin - 1;
%         p_switch_bin(bin_i+1) = sum(p_switch(devt_i_start:devt_i_end) .* n_delta_events(devt_i_start:devt_i_end)) / sum(n_delta_events(devt_i_start:devt_i_end));
%         % p_switch_bin(bin_i+1) = nanmean(p_switch(devt_i_start:devt_i_end));
%     end
%     p_switch = p_switch_bin;
% end

end
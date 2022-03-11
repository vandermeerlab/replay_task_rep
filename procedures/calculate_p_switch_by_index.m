function [p_switch] = calculate_p_switch_by_index(cfg_in, SWR_data)

% Calculate the probablity of switching across deltas of event index

cfg_def = [];

cfg = ProcessConfig(cfg_def,cfg_in);

L2_L1_count = 0;
R2_R1_count = 0;
R2_L1_count = 0;
L2_R1_count = 0;

for evt_i = 1:length(SWR_data)-1
    if ~isnan(SWR_data(evt_i+1)) && ~isnan(SWR_data(evt_i))
        if SWR_data(evt_i+1) > 0 && SWR_data(evt_i) > 0
            L2_L1_count = L2_L1_count + 1;
        elseif SWR_data(evt_i+1) < 0 && SWR_data(evt_i) < 0
            R2_R1_count = R2_R1_count + 1;
        elseif SWR_data(evt_i+1) < 0 && SWR_data(evt_i) > 0
            R2_L1_count = R2_L1_count + 1;
        elseif SWR_data(evt_i+1) > 0 && SWR_data(evt_i) < 0
            L2_R1_count = L2_R1_count + 1;
        end
    end
end
L1_count = L2_L1_count + R2_L1_count;
R1_count = L2_R1_count + R2_R1_count;

if L1_count == 0
    P_R2_L1 = 0;
    P_L2_L1 = 0;
else
    P_R2_L1 = R2_L1_count / L1_count;
    P_L2_L1 = L2_L1_count / L1_count;
end

if R1_count == 0
    P_L2_R1 = 0;
    P_R2_R1 = 0;
else
    P_L2_R1 = L2_R1_count / R1_count;
    P_R2_R1 = R2_R1_count / R1_count;
end

p_switch = (P_R2_L1 + P_L2_R1) - (P_L2_L1 + P_R2_R1);

end
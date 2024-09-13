%%
first_half_nums = zeros(1, length(out));
second_half_nums = zeros(1, length(out));

alt_oppo_count = zeros(1, length(out));
L_oppo_count = zeros(1, length(out));
R_oppo_count = zeros(1, length(out));

alt_total_count = zeros(1, length(out));
L_total_count = zeros(1, length(out));
R_total_count = zeros(1, length(out));

alt_shuf_z = zeros(1, length(out));
L_shuf_z = zeros(1, length(out));
R_shuf_z = zeros(1, length(out));

oppo_replay_bias = zeros(1, length(out));

sig_cutoff = 0.05;

for s_i = 1:length(out)
    switch_idx = out{s_i}.switch_idx;
    first_half_nums(s_i) = length(out{s_i}.behav_sequence(1:switch_idx));
    second_half_nums(s_i) = length(out{s_i}.behav_sequence(switch_idx+1:end));
    oppo_bias = [];
    for c_i = 1:length(out{s_i}.contigency)
        shuf_z = [];
        if c_i == 1
            block_num = length(out{s_i}.behav_sequence(1:switch_idx));
        else
            block_num = length(out{s_i}.behav_sequence(switch_idx+1:end));
        end
        % Only use half-sessions with > 10 trials
        if block_num > 10
            for t_i = 1:length(out{s_i}.trial_iv.tend)
                if t_i <= length(out{s_i}.trial_iv.tend)-1
                    trial_SWRs = find(out{s_i}.tvec > out{s_i}.trial_iv.tend(t_i) & ...
                        out{s_i}.tvec < out{s_i}.trial_iv.tstart(t_i+1));
                else
                    trial_SWRs = find(out{s_i}.tvec > out{s_i}.trial_iv.tend(t_i) & ...
                        out{s_i}.tvec < out{s_i}.TimeOffTrack);
                end
                for re_i = 1:length(trial_SWRs)
                    SWR_index = trial_SWRs(re_i);
                   if out{s_i}.shuf_perc_odds(SWR_index) < sig_cutoff || ...
                           out{s_i}.shuf_perc_odds(SWR_index) > 1-sig_cutoff

                        if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                            alt_total_count(s_i) = alt_total_count(s_i)+1;
                        elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
                            L_total_count(s_i) = L_total_count(s_i)+1;
                        else
                            R_total_count(s_i) = R_total_count(s_i)+1;
                        end

                        if strcmp(out{s_i}.contigency{c_i}, 'Left')
                            oppo_bias = [oppo_bias, -out{s_i}.shuf_z_diff(SWR_index)];
                        elseif strcmp(out{s_i}.contigency{c_i}, 'Right')
                            oppo_bias = [oppo_bias, out{s_i}.shuf_z_diff(SWR_index)];
                        end

                        shuf_z = [shuf_z, out{s_i}.shuf_z_diff(SWR_index)];

%                         if out{s_i}.actual_diff(SWR_index) > 0
                        if out{s_i}.shuf_perc_odds(SWR_index) > 1-sig_cutoff
                            % Left is more represented in this SWR
                            if strcmp(out{s_i}.behav_sequence{t_i}, 'R')
                                if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                                    alt_oppo_count(s_i) = alt_oppo_count(s_i)+1;
                                elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
                                    L_oppo_count(s_i) = L_oppo_count(s_i)+1;
                                else
                                    R_oppo_count(s_i) = R_oppo_count(s_i)+1;
                                end
                            end
                        else
                            % Right is more represented in this SWR
                            if strcmp(out{s_i}.behav_sequence{t_i}, 'L')
                                if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                                    alt_oppo_count(s_i) = alt_oppo_count(s_i)+1;
                                elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
                                    L_oppo_count(s_i) = L_oppo_count(s_i)+1;
                                else
                                    R_oppo_count(s_i) = R_oppo_count(s_i)+1;
                                end
                            end
                        end
                    end
                end
            end
        end
        if strcmp(out{s_i}.contigency{c_i}, 'Alt')
            alt_shuf_z(s_i) = nanmedian(shuf_z);
        elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
            L_shuf_z(s_i) = nanmedian(shuf_z);
        else
            R_shuf_z(s_i) = nanmedian(shuf_z);
        end
        oppo_replay_bias(s_i) = nanmedian(oppo_bias);
    end
end    

%%
alt_oppo_prop = NaN(1, length(out));
L_oppo_prop = NaN(1, length(out));
R_oppo_prop = NaN(1, length(out));

for s_i = 1:length(out)
    if alt_total_count(s_i) ~= 0
        alt_oppo_prop(s_i) = alt_oppo_count(s_i) / alt_total_count(s_i);
    end
    if L_total_count(s_i) ~= 0
        L_oppo_prop(s_i) = L_oppo_count(s_i) / L_total_count(s_i);
    end
    if R_total_count(s_i) ~= 0
        R_oppo_prop(s_i) = R_oppo_count(s_i) / R_total_count(s_i);
    end
end

%%
types = {'Alt (sig.)', 'Left (sig.)', 'Right (sig.)'};
% oppo_props = {alt_oppo_prop, L_oppo_prop, R_oppo_prop};
oppo_props = {alt_shuf_z, L_shuf_z, R_shuf_z};
mean_oppo_prop = zeros(length(oppo_props), 1);
sem_oppo_prop = zeros(length(oppo_props), 1);

for d_i = 1:length(oppo_props)
    mean_oppo_prop(d_i) = nanmedian(oppo_props{d_i});
    sem_oppo_prop(d_i) = nanstd(oppo_props{d_i}) / sqrt(sum(~isnan(oppo_props{d_i})));
end

ylim = [-1, 1];
dx = 0.05;
dy = 0.1;
fs = 12;
% y_label = 'proportion of opposite-side replay';
y_label = 'z-scored replay bias';

x = dx * (1:length(oppo_props));
xpad = 0.05;
h = errorbar(x, mean_oppo_prop, sem_oppo_prop, 'LineStyle', 'none', 'LineWidth', 2);
set(h, 'Color', 'k');
hold on;
plot(x, mean_oppo_prop, '.k', 'MarkerSize', 20);
set(gca, 'XTick', x, 'YTick', [ylim(1):dy:ylim(2)], 'XTickLabel', types, ...
    'XLim', [x(1)-xpad x(end)+xpad], 'YLim', [ylim(1), ylim(2)], 'FontSize', fs, ...
    'LineWidth', 1, 'TickDir', 'out');
box off;
plot([x(1)-xpad x(end)+xpad], [0 0], '--k', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
ylabel(y_label);
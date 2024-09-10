%%
first_half_nums = zeros(1, length(out));
second_half_nums = zeros(1, length(out));

alt_oppo_count = zeros(1, length(out));
repeat_oppo_count = zeros(1, length(out));

alt_total_count = zeros(1, length(out));
repeat_total_count = zeros(1, length(out));

sig_cutoff = 0.05;

for s_i = 1:length(out)
    switch_idx = out{s_i}.switch_idx;
    first_half_nums(s_i) = length(out{s_i}.behav_sequence(1:switch_idx));
    second_half_nums(s_i) = length(out{s_i}.behav_sequence(switch_idx+1:end));
    for c_i = 1:length(out{s_i}.contigency)
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
                   if out{s_i}.shuf_perc_diff(SWR_index) < sig_cutoff || ...
                           out{s_i}.shuf_perc_diff(SWR_index) > 1-sig_cutoff

                        if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                            alt_total_count(s_i) = alt_total_count(s_i)+1;
                        else
                            repeat_total_count(s_i) = repeat_total_count(s_i)+1;
                        end

                        if out{s_i}.actual_diff(SWR_index) > 0
                            % Left is more represented in this SWR
                            if strcmp(out{s_i}.behav_sequence{t_i}, 'R')
                                if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                                    alt_oppo_count(s_i) = alt_oppo_count(s_i)+1;
                                else
                                    repeat_oppo_count(s_i) = repeat_oppo_count(s_i)+1;
                                end
                            end
                        else
                            % Right is more represented in this SWR
                            if strcmp(out{s_i}.behav_sequence{t_i}, 'L')
                                if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                                    alt_oppo_count(s_i) = alt_oppo_count(s_i)+1;
                                else
                                    repeat_oppo_count(s_i) = repeat_oppo_count(s_i)+1;
                                end
                            end
                        end
                   end
                end
            end
        end
    end
end    

%%
alt_oppo_prop = NaN(1, length(out));
repeat_oppo_prop = NaN(1, length(out));

for s_i = 1:length(out)
    if alt_total_count(s_i) ~= 0
        alt_oppo_prop(s_i) = alt_oppo_count(s_i) / alt_total_count(s_i);
    end
    if repeat_total_count(s_i) ~= 0
        repeat_oppo_prop(s_i) = repeat_oppo_count(s_i) / repeat_total_count(s_i);
    end
end

%%
types = {'Alt (sig.)', 'Repeat (sig.)'};
oppo_props = {alt_oppo_prop, repeat_oppo_prop};
mean_oppo_prop = zeros(length(oppo_props), 1);
sem_oppo_prop = zeros(length(oppo_props), 1);

for d_i = 1:length(oppo_props)
    mean_oppo_prop(d_i) = nanmean(oppo_props{d_i});
    sem_oppo_prop(d_i) = nanstd(oppo_props{d_i}) / sqrt(sum(~isnan(oppo_props{d_i})));
end

ylim = [0.3, 0.7];
dx = 0.05;
dy = 0.1;
fs = 12;
y_label = 'proportion of opposite-side replay';

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
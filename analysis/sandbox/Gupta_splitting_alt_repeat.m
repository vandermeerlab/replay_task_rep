%%

alt_total_count = NaN(1, length(out));
L_total_count = NaN(1, length(out));
R_total_count = NaN(1, length(out));

for s_i = 1:length(out)
    switch_idx = out{s_i}.switch_idx;
    for c_i = 1:length(out{s_i}.contigency)
        if c_i == 1
            block_num = length(out{s_i}.behav_sequence(1:switch_idx));
        else
            block_num = length(out{s_i}.behav_sequence(switch_idx+1:end));
        end
        % Only use half-sessions with > 10 trials
        if block_num > 10
            if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                alt_total_count(s_i) = block_num;
            elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
                L_total_count(s_i) = block_num;
            else
                R_total_count(s_i) = block_num;
            end
        end
    end
end

%%
alt_errs_count = NaN(1, length(out));
L_errs_count = NaN(1, length(out));
R_errs_count = NaN(1, length(out));

for s_i = 1:length(out)
    alt_errs_count(s_i) = alt_total_count(s_i) - alt_correct_count(s_i);
    L_errs_count(s_i) = L_total_count(s_i) - L_correct_count(s_i);
    R_errs_count(s_i) = R_total_count(s_i) - R_correct_count(s_i);
end

repeat_err_count = [L_errs_count, R_errs_count];
repeat_err_count(repeat_err_count == 0) = NaN;

%%
alt_splitting = NaN(1, length(out));
L_splitting = NaN(1, length(out));
R_splitting = NaN(1, length(out));

for s_i = 1:length(out)
    switch_idx = out{s_i}.switch_idx;
    for c_i = 1:length(out{s_i}.contigency)
        FR_diff = mean(out{s_i}.FR_diff{c_i}, 'omitnan');
        if c_i == 1
            block_num = length(out{s_i}.behav_sequence(1:switch_idx));
        else
            block_num = length(out{s_i}.behav_sequence(switch_idx+1:end));
        end
        % Only use half-sessions with > 10 trials
        if block_num > 10
            if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                alt_total_count(s_i) = block_num;
                alt_splitting(s_i) = FR_diff;
            elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
                L_total_count(s_i) = block_num;
                if L_errs_count(s_i) > 2
                    L_splitting(s_i) = FR_diff;
                end
            else
                R_total_count(s_i) = block_num;
                if R_errs_count(s_i) > 2
                    R_splitting(s_i) = FR_diff;
                end
            end
        end
    end
end

%%
%%
types = {'Alt', 'R or L'};
single_splitting = {alt_splitting, [L_splitting, R_splitting]};
mean_splittng = zeros(length(single_splitting), 1);
sem_splittng = zeros(length(single_splitting), 1);

for d_i = 1:length(single_splitting)
    mean_splittng(d_i) = nanmean(single_splitting{d_i});
    sem_splittng(d_i) = nanstd(single_splitting{d_i}) / sqrt(sum(~isnan(single_splitting{d_i})));
end

ylim = [0.4, 1.1];
dx = 0.05;
dy = 0.2;
fs = 24;
y_label = '$\left|FR_{left} - FR_{right}\right|$';

x = dx * (1:length(single_splitting));
xpad = 0.05;
h = errorbar(x, mean_splittng, sem_splittng, 'LineStyle', 'none', 'LineWidth', 2);
set(h, 'Color', 'k');
hold on;
plot(x, mean_splittng, '.k', 'MarkerSize', 20);
set(gca, 'XTick', x, 'YTick', [ylim(1):dy:ylim(2)], 'XTickLabel', types, ...
    'XLim', [x(1)-xpad x(end)+xpad], 'YLim', [ylim(1), ylim(2)], 'FontSize', fs, ...
    'LineWidth', 1, 'TickDir', 'out');
box off;
plot([x(1)-xpad x(end)+xpad], [0 0], '--k', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
ylabel(y_label, 'Interpreter','latex');
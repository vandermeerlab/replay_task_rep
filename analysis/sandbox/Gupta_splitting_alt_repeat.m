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
                alt_splitting(s_i) = FR_diff;
            elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
                if L_errs_count(s_i) > 2
                    L_splitting(s_i) = FR_diff;
                end
            else
                if R_errs_count(s_i) > 2
                    R_splitting(s_i) = FR_diff;
                end
            end
        end
    end
end

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

%%
alt_ensemble_splitting = NaN(1, length(out));
L_ensemble_splitting = NaN(1, length(out));
R_ensemble_splitting = NaN(1, length(out));

for s_i = 1:length(out)
    switch_idx = out{s_i}.switch_idx;
    for c_i = 1:length(out{s_i}.contigency)
        avg_preCP_correct = max(out{s_i}.left_preCP_correct{c_i}, out{s_i}.right_preCP_correct{c_i});
        if c_i == 1
            block_num = length(out{s_i}.behav_sequence(1:switch_idx));
        else
            block_num = length(out{s_i}.behav_sequence(switch_idx+1:end));
        end
        % Only use half-sessions with > 10 trials
        if block_num > 10
            if strcmp(out{s_i}.contigency{c_i}, 'Alt')
                alt_ensemble_splitting(s_i) = avg_preCP_correct;
            elseif strcmp(out{s_i}.contigency{c_i}, 'Left')
                if L_errs_count(s_i) > 2
                    L_ensemble_splitting(s_i) = avg_preCP_correct;
                end
            else
                if R_errs_count(s_i) > 2
                    R_ensemble_splitting(s_i) = avg_preCP_correct;
                end
            end
        end
    end
end

%%
types = {'Alt', 'R or L'};
ensemble_splitting = {alt_ensemble_splitting, [L_ensemble_splitting, R_ensemble_splitting]};
mean_splittng = zeros(length(ensemble_splitting), 1);
sem_splittng = zeros(length(ensemble_splitting), 1);

for d_i = 1:length(ensemble_splitting)
    mean_splittng(d_i) = nanmean(ensemble_splitting{d_i});
    sem_splittng(d_i) = nanstd(ensemble_splitting{d_i}) / sqrt(sum(~isnan(ensemble_splitting{d_i})));
end

ylim = [0.5, 0.8];
dx = 0.05;
dy = 0.05;
fs = 24;
y_label = 'ensemble splitter strength';

x = dx * (1:length(ensemble_splitting));
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

%%
x = single_splitting{1};
y = ensemble_splitting{1};

x = x(~isnan(x));
y = y(~isnan(y));

figure;
p = polyfit(x, y, 1);
predicted_y = polyval(p, x);

[corr, p_val] = corrcoef(x, y);
scatter(x, y, 'filled', 'MarkerFaceColor', [105/255 105/255 105/255]);
hold on;
plot(x, predicted_y, '--', 'Color', [198/255 113/255 113/255], 'LineWidth', 2); 

xlabel('single-cell splitter strength')
ylabel('ensemble splitter strength')
title(['r = ', num2str(corr(2, 1)), '; p-value = ', num2str(p_val(2, 1))])
% title(['r = ', num2str(bias_corr(2, 1))])

set(gca, 'FontSize', 24)
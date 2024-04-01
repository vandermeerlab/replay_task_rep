food_idx = ~logical(data.all.all.this_type);
water_idx = logical(data.all.all.this_type);

replay_bias_scores = zeros(length(data.all.all.median_z), 1);
replay_bias_scores(food_idx) = -data.all.all.median_z(food_idx);
replay_bias_scores(water_idx) = data.all.all.median_z(water_idx);

avg_preCP_correct = (left_preCP_correct + right_preCP_correct) / 2;

n_neurons = data.all.all.n_neurons;
FR_diff = data.all.all.FR_diff;
FR_diff_zscore = data.all.all.FR_diff_zscore;

%%
x = avg_preCP_correct;
y = replay_bias_scores;

figure;
p = polyfit(x, y, 1);
predicted_bias = polyval(p, x);

[bias_corr, p_val] = corrcoef(x, y);
scatter(x, y, 'filled', 'MarkerFaceColor', [105/255 105/255 105/255]);
hold on;
plot(x, predicted_bias, '--', 'Color', [198/255 113/255 113/255], 'LineWidth', 2); 

xlabel('splitter strength')
ylabel('replay bias')
% title(['r = ', num2str(bias_corr(2, 1)), '; p-value = ', num2str(p_val(2, 1))])
title(['r = ', num2str(bias_corr(2, 1))])
set(gca, 'FontSize', 24)

%%
x = FR_diff;
y = replay_bias_scores;

figure;
p = polyfit(x, y, 1);
predicted_bias = polyval(p, x);

[bias_corr, p_val] = corrcoef(x, y);
scatter(x, y, 'filled', 'MarkerFaceColor', [105/255 105/255 105/255]);
hold on;
plot(x, predicted_bias, '--', 'Color', [198/255 113/255 113/255], 'LineWidth', 2); 

xlabel('$\left|FR_{left} - FR_{right}\right|$ ', 'Interpreter','latex')
ylabel('replay bias')
% title(['r = ', num2str(bias_corr(2, 1)), '; p-value = ', num2str(p_val(2, 1))])
title(['r = ', num2str(bias_corr(2, 1))])
set(gca, 'FontSize', 24)

%%
avg_correct_alls = {avg_preCP_correct, avg_postCP_correct, avg_entire_correct};
accuracy_titles = {'pre-CP', 'post-CP', 'entire'};
for a_i = 1:length(avg_correct_alls)
    figure;
    avg_correct = avg_correct_alls{a_i};
    scatter(avg_correct(food_idx), -data.all.all.median_z(food_idx), 'filled', 'red');
    hold on;
    scatter(avg_correct(water_idx), data.all.all.median_z(water_idx), 'filled', 'blue');
    xlabel('Decoding accuracy')
    ylabel('opposite-to-behavior bias in replay')
    legend({'food-restricted', 'water-restricted'})
%     xlim([0.4, 1])
    title(accuracy_titles{a_i})
    set(gca, 'FontSize', 14)
end

%% replay bias per subject
subj_ids = {'R042', 'R044', 'R50', 'R064'};
subj_data_idx = {1:5, 6:7, 8:13, 14:19};
subj_sess_idx = {2:6, 4:5, 1:6, 1:6};
subj_colors = {[143/255 188/255 143/255], [176/255 196/255 222/255], ...
    [216/255 191/255 216/255], [238/255 213/255 210/255]};

figure;
set(gcf, 'Position', [734 355 615 511]);

x = 1:6;
for s_i = 1:length(subj_ids)
    y = NaN(1, 6);
    y(subj_sess_idx{s_i}) = replay_bias_scores(subj_data_idx{s_i});
    plot(x, y, '.-', 'Color', subj_colors{s_i}, 'DisplayName',subj_ids{s_i}, ...
        'MarkerSize', 20, 'LineWidth', 2);
    hold on;
end
xlabel('day');
ylabel('replay bias');
xpad = 0.25;
xlim = [x(1)-xpad x(end)+xpad];

ypad = 0.6;
ylim = [-0.9, 1.5];
yt = ylim(1):ypad:ylim(2);
ytl = {ylim(1), '', (ylim(1) + ylim(2)) / 2, '', ylim(2)};

set(gca, 'XTick', x, 'YTick', yt, 'YTickLabel', ytl, ...
    'XTickLabel', x, 'XLim', xlim, 'YLim', ylim, ...
    'FontSize', 24,'LineWidth', 1, 'TickDir', 'out');
legend('Location','southwest');
box off;

%% splitter strength per subject
figure;
set(gcf, 'Position', [734 355 615 511]);

x = 1:6;
for s_i = 1:length(subj_ids)
    y = NaN(1, 6);
    y(subj_sess_idx{s_i}) = avg_preCP_correct(subj_data_idx{s_i});
    plot(x, y, '.-', 'Color', subj_colors{s_i}, 'DisplayName',subj_ids{s_i}, ...
        'MarkerSize', 20, 'LineWidth', 2);
    hold on;
end
xlabel('day');
ylabel('splitter strength');
xpad = 0.25;
xlim = [x(1)-xpad x(end)+xpad];

ypad = 0.15;
ylim = [0.4, 1.0];
yt = ylim(1):ypad:ylim(2);
ytl = {ylim(1), '', (ylim(1) + ylim(2)) / 2, '', ylim(2)};

set(gca, 'XTick', x, 'YTick', yt, 'YTickLabel', ytl, ...
    'XTickLabel', x, 'XLim', xlim, 'YLim', ylim, ...
    'FontSize', 24,'LineWidth', 1, 'TickDir', 'out');
legend('Location','southwest');
box off;

%% n_neurons per subject
figure;
set(gcf, 'Position', [734 355 615 511]);

x = 1:6;
for s_i = 1:length(subj_ids)
    y = NaN(1, 6);
    y(subj_sess_idx{s_i}) = n_neurons(subj_data_idx{s_i});
    plot(x, y, '.-', 'Color', subj_colors{s_i}, 'DisplayName',subj_ids{s_i}, ...
        'MarkerSize', 20, 'LineWidth', 2);
    hold on;
end
xlabel('day');
ylabel('number of neurons');
xpad = 0.25;
xlim = [x(1)-xpad x(end)+xpad];

ypad = 15;
ylim = [20, 80];
yt = ylim(1):ypad:ylim(2);
ytl = {ylim(1), '', (ylim(1) + ylim(2)) / 2, '', ylim(2)};

set(gca, 'XTick', x, 'YTick', yt, 'YTickLabel', ytl, ...
    'XTickLabel', x, 'XLim', xlim, 'YLim', ylim, ...
    'FontSize', 24,'LineWidth', 1, 'TickDir', 'out');
legend('Location','southwest');
box off;
%% Plot odd ratios against real world time
SWR_data = actual_L_R_diff;
SWR_data_L = NaN(size(SWR_data));
SWR_data_L(SWR_data > 0) = SWR_data(SWR_data > 0);
SWR_data_R = NaN(size(SWR_data));
SWR_data_R(SWR_data < 0) = SWR_data(SWR_data < 0);

stem(Q_SWR.tvec, SWR_data_L, 'Color', 'r');
hold on;
stem(Q_SWR.tvec, SWR_data_R, 'Color', [0 0.4470 0.7410]);
xline(ExpKeys.prerecord(2), '--k', 'Pre-task end');
% xline(metadata.taskvars.rest_iv.tend(end), '--r', 'Task end');
xline(ExpKeys.postrecord(1), '--k', 'Post-task start');

xlim([Q_SWR.tvec(1) - 100, Q_SWR.tvec(end) + 100]);
ylim([-1.25, 1.25]);
xlabel('Time (s)')
ylabel('pL - pR')
set(gca,'FontSize', 18)

%% Plot autocorrelation of pL - pR against SWR event indices
n_lags = 50;
SWR_data = actual_L_R_diff;

[acf,lags] = xcorr(SWR_data(~isnan(SWR_data) & ~isinf(SWR_data)), n_lags, 'coeff');
plot(lags, acf);
hold on;
ylabel('pL - pR');
xlabel('lags between SWR indices');
set(gca,'FontSize', 18)

%% Shuffling SWR indices
n_shuffles = 1000;
acf_shuffles = zeros(n_shuffles, 2*n_lags+1);

for s_i = 1:n_shuffles
    shuffle_indices = randperm(length(pL));
    s_SWR_data = SWR_data(shuffle_indices);
    [acf, lags] = xcorr(s_SWR_data(~isnan(s_SWR_data) & ~isinf(s_SWR_data)), n_lags, 'coeff');
    acf_shuffles(s_i, :) = acf;
end

shuffle_mean = mean(acf_shuffles, 1);
u_bound = prctile(acf_shuffles, 97.5, 1) - shuffle_mean;
l_bound =shuffle_mean - prctile(acf_shuffles, 2.5, 1);
h = shadedErrorBar(lags, shuffle_mean, [u_bound; l_bound]);

%% Select significant events
L_sig_odds = out.shuf_perc_odds >= 0.95;
R_sig_odds = out.shuf_perc_odds <= 0.05;

L_sig_diff = out.shuf_perc_diff >= 0.95;
R_sig_diff = out.shuf_perc_diff <= 0.05;

%%
Q_SWR = MakeQfromS(cfg_Q,expComb.decS);
Q_SWR.data = Q_SWR.data(:,1:2:end);
Q_SWR.tvec = Q_SWR.tvec(1:2:end);

% need to do some selection on minimum number of cells?
nActiveCells = sum(Q_SWR.data > 0);
keep = nActiveCells >= cfg.nMinNeurons;
Q_SWR.tvec = Q_SWR.tvec(keep); Q_SWR.data = Q_SWR.data(:,keep);

%% Plot L as red and R as blue against real world time
SWR_data = actual_L_R_diff;
L_or_R_sig = L_sig | R_sig;
sig_time = Q_SWR.tvec(L_or_R_sig);
L_sig_discrete = ones(size(SWR_data(L_sig)));
R_sig_discrete = ones(size(SWR_data(R_sig)));

stem(Q_SWR.tvec(L_sig), L_sig_discrete, 'Marker', 'none', 'Color', 'r');
hold on;
stem(Q_SWR.tvec(R_sig), R_sig_discrete, 'Marker', 'none', 'Color', [0 0.4470 0.7410]);

xline(ExpKeys.prerecord(2), '--k', 'Pre-task end', 'HandleVisibility','off');
% xline(metadata.taskvars.rest_iv.tend(end), '--r', 'Task end');
xline(ExpKeys.postrecord(1), '--k', 'Post-task start', 'HandleVisibility','off');
legend('Left','Right', '', '');

xlim([sig_time(1) - 100, sig_time(end) + 100]);
ylim([-0.25, 1.25]);
xlabel('Time (s)')
set(gca, 'ytick',[], 'FontSize', 18)

%% Plot autocorrelation of odd ratos against (significant) SWR event indices (pre-task, ITI, post-task)
bev_state_titles = {'Pre-task', 'ITI', 'Post-task'};
L_or_R_sig = L_sig | R_sig;
SWR_data = actual_odds;

for i = 1:size(bev_state_titles, 2)
    bev_state_title = bev_state_titles{i};
    if i == 1
        [~,keep] = restrict_idx(Q_SWR,0,ExpKeys.prerecord(2));
    elseif i == 2
        [~,keep] = restrict_idx(Q_SWR,metadata.taskvars.rest_iv);
    else
        [~,keep] = restrict_idx(Q_SWR,ExpKeys.postrecord(1),Inf);
    end
    SWR_data_state = SWR_data(keep);
    [acf,lags] = xcorr(SWR_data_state(~isnan(SWR_data_state) & ~isinf(SWR_data_state)), 10, 'coeff');
    subplot(1, 3, i);
    plot(lags, acf);
    ylabel('log odd ratio');
    xlabel('lags between SWR indices');
    title(bev_state_title);
    set(gca,'FontSize', 18)
end

%% Plot the distribution Inter-SWRs-interval across sessions
SWR_time_diffs = [];
for sess_i = 1:length(out)
    SWR_times = out{sess_i}.tvec;
    SWR_time_diffs = [SWR_time_diffs; diff(SWR_times)];
end

figure;
histogram(log10(SWR_time_diffs));
xlabel('log10 (time elapsed since last SWR)')
ylabel('Proportion')
set(gca, 'ytick',[], 'FontSize', 18)
title('Distribution of Inter-SWRs-interval')

%% Plot Inter-SWRs-interval as a function of time.
sig_events_only = 0;
bin_size = 0.1;
bin_egdes = -1:bin_size:2;
devt_max = 5;

p_switch = cell(1, length(out));
p_switch_baseline = zeros(1, length(out));
adj_p_switch = [];
num_events_bin_mat = [];

for sess_i = 1:length(out)
    SWR_data = out{sess_i}.actual_pL - out{sess_i}.actual_pR;
    SWR_times = out{sess_i}.tvec;

    if sig_events_only
        L_sig_odds = out{sess_i}.shuf_perc_odds >= 0.95;
        R_sig_odds = out{sess_i}.shuf_perc_odds <= 0.05;
        SWR_data = SWR_data(L_sig_odds | R_sig_odds);
        SWR_times = SWR_times(L_sig_odds | R_sig_odds);
    end
    SWR_data = SWR_data(~isnan(SWR_data));
    SWR_times = SWR_times(~isnan(SWR_data));

    cfg_sw = [];
    cfg_sw.bin_egdes = bin_egdes;
    cfg_sw.devt_max = devt_max;
    [p_switch{sess_i}, num_events_bin{sess_i}] = calculate_p_switch_by_time(cfg_sw, SWR_data, SWR_times);
    num_events_bin_mat = [num_events_bin_mat; num_events_bin{sess_i}];

    % Calculate baseline: 2*pL*pR
    pL = sum(SWR_data > 0) / length(SWR_data);
    p_switch_baseline(sess_i) = 2 * pL * (1-pL);

    adj_p_switch = [adj_p_switch; p_switch{sess_i} - p_switch_baseline(sess_i)];
end

%% Shuffling SWR indices
n_shuffles = 1000;
p_switch_shuffles = cell(1, length(out));
adj_p_switch_shuffles = cell(1, length(out));

for sess_i = 1:length(out)

    SWR_data = out{sess_i}.actual_pL - out{sess_i}.actual_pR;
    SWR_times = out{sess_i}.tvec;

    if sig_events_only
        L_sig_odds = out{sess_i}.shuf_perc_odds >= 0.95;
        R_sig_odds = out{sess_i}.shuf_perc_odds <= 0.05;
        SWR_data = SWR_data(L_sig_odds | R_sig_odds);
        SWR_times = SWR_times(L_sig_odds | R_sig_odds);
    end
    SWR_data = SWR_data(~isnan(SWR_data));
    SWR_times = SWR_times(~isnan(SWR_data));

    cfg_sw = [];
    cfg_sw.bin_egdes = bin_egdes;
    cfg_sw.devt_max = devt_max;

    for s_i = 1:n_shuffles
        shuffle_indices = randperm(length(SWR_data));
        s_SWR_data = SWR_data(shuffle_indices);
        [s_p_switch] = calculate_p_switch_by_time(cfg_sw, s_SWR_data, SWR_times);

        p_switch_shuffles{sess_i}(s_i, :) = s_p_switch;
        adj_p_switch_shuffles{sess_i}(s_i, :) = s_p_switch - p_switch_baseline(sess_i);
    end
end

%% Plotting probablity of switching with baseline for each session
for sess_i = 1:length(out)
    t_diffs_x = bin_egdes(1:end-1) + bin_size / 2;

    figure;
    plot(t_diffs_x, p_switch{sess_i}, '.-r'); hold on;
    yline(p_switch_baseline(sess_i), '--k', 'Baseline'); hold on;

    u_bound = prctile(p_switch_shuffles{sess_i}, 97.5, 1) - mean(p_switch_shuffles{sess_i}, 1);
    l_bound = mean(p_switch_shuffles{sess_i}, 1) - prctile(p_switch_shuffles{sess_i}, 2.5, 1);
    h = shadedErrorBar(t_diffs_x, mean(p_switch_shuffles{sess_i}, 1), [u_bound;l_bound]);

    xlim([-1, 2]);
    ylim([-0.25, 1.25]);
    xlabel('log10 (time elapsed since last SWR)')
    ylabel('P(interleave)')
    set(gca,'FontSize', 18)

    if sig_events_only
        title(sprintf('Significant events by odd ratio Session %d', sess_i));
    else
        title(sprintf('All events Session %d', sess_i));
    end
    % saveas(gcf, sprintf('all_events_%d.jpg', sess_i));
end

%%
adj_p_switch_shuffles_mat = zeros(1000, length(cfg_sw.bin_egdes)-1);
for i = 1:1000
    shuffles_mat = [];
    for sess_i = 1:length(out)
        if isempty(shuffles_mat)
            shuffles_mat = adj_p_switch_shuffles{sess_i}(i, :);
        else
            shuffles_mat = [shuffles_mat; adj_p_switch_shuffles{sess_i}(i, :)];
        end
    end
    adj_p_switch_shuffles_mat(i, :) = nanmean(shuffles_mat, 1);
end

%% Plotting adjusted (by independent baseline) probablity of switching across sessions
adj_p_switch_m = nanmean(adj_p_switch, 1);
adj_p_switch_sem = nanstd(adj_p_switch, 1) / sqrt(size(adj_p_switch, 1));

x = bin_egdes(1:end-1) + bin_size / 2;
xpad = 0.05;
h = errorbar(x, adj_p_switch_m, adj_p_switch_sem, 'LineWidth', 1.5); hold on;
set(h, 'Color', 'k');

adj_p_switch_shuffles_m = mean(adj_p_switch_shuffles_mat, 1);
u_bound = prctile(adj_p_switch_shuffles_mat, 97.5, 1) - adj_p_switch_shuffles_m;
l_bound = adj_p_switch_shuffles_m - prctile(adj_p_switch_shuffles_mat, 2.5, 1);
sh = shadedErrorBar(x, adj_p_switch_shuffles_m, [u_bound;l_bound]);

hold on;
plot(x, adj_p_switch_m, '.k', 'MarkerSize', 20);
set(gca, 'XLim', [x(1)-xpad x(end)+xpad], 'YLim', [-0.25, 0.25], 'FontSize', 18, ...
    'LineWidth', 1, 'TickDir', 'out');
box off;
plot([x(1)-xpad x(end)+xpad], [0 0], '--k', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);

xlabel('log10 (time elapsed since last event)')
ylabel('adjusted P(interleave)')

if sig_events_only
    title('Significant events by odd ratio delta 1');
else
    title('All events delta 5');
end
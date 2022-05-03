%% Calculate adjusted (by independent baseline) probablity of switching and baseline by delta event index
devt_max = 100;
% n_bins = devt_max / devt_nbin;
sig_events_only = 1;

p_switch = cell(1, length(out));
p_switch_baseline = zeros(1, length(out));
adj_p_switch = zeros(length(out), devt_max);

for sess_i = 1:length(out)
    SWR_data = out{sess_i}.actual_pL - out{sess_i}.actual_pR;

    if sig_events_only
        L_sig_odds = out{sess_i}.shuf_perc_odds >= 0.95;
        R_sig_odds = out{sess_i}.shuf_perc_odds <= 0.05;
        SWR_data = SWR_data(L_sig_odds | R_sig_odds);
    end
    SWR_data = SWR_data(~isnan(SWR_data));

    cfg_switch = [];
    cfg_switch.devt_max = devt_max;
%     cfg_switch.devt_nbin = devt_nbin;
    [p_switch{sess_i}] = calculate_p_switch_by_index(cfg_switch, SWR_data);
    % Calculate baseline: 2*pL*pR
    pL = sum(SWR_data > 0) / length(SWR_data);
    p_switch_baseline(sess_i) = 2 * pL * (1-pL);

    adj_p_switch(sess_i, :) = p_switch{sess_i} - p_switch_baseline(sess_i);
end

%% Shuffling SWR indices
n_shuffles = 1000;
p_switch_shuffles = cell(1, length(out));
adj_p_switch_shuffles = cell(1, length(out));

for sess_i = 1:length(out)
    cfg_switch = [];
    cfg_switch.devt_max = devt_max;
%     cfg_switch.devt_nbin = devt_nbin;

    p_switch_shuffles{sess_i} = zeros(n_shuffles, devt_max);
    SWR_data = out{sess_i}.actual_pL - out{sess_i}.actual_pR;

    if sig_events_only
        L_sig_odds = out{sess_i}.shuf_perc_odds >= 0.95;
        R_sig_odds = out{sess_i}.shuf_perc_odds <= 0.05;
        SWR_data = SWR_data(L_sig_odds | R_sig_odds);
    end
    SWR_data = SWR_data(~isnan(SWR_data));

    for s_i = 1:n_shuffles
        shuffle_indices = randperm(length(SWR_data));
        s_SWR_data = SWR_data(shuffle_indices);

        % SWR times for significant events should stay constant
        s_p_switch = calculate_p_switch_by_index(cfg_switch, s_SWR_data);

        p_switch_shuffles{sess_i}(s_i, :) = s_p_switch;
        adj_p_switch_shuffles{sess_i}(s_i, :) = s_p_switch - p_switch_baseline(sess_i);
    end
end

%% Calculate adjusted probablity of switching for behavioral sequence (L or R)
devt_max_behav = 10;

p_switch_behav = cell(1, length(out));
p_switch_baseline_behav = zeros(1, length(out));
adj_p_switch_behav = zeros(length(out), devt_max_behav);

for sess_i = 1:length(out)
    behav_sequence = out{sess_i}.behav_sequence;
    behav_data = [];
    for b_i = 1:length(behav_sequence)
        if strcmp(behav_sequence{b_i}, 'L')
            behav_data = [behav_data, 1];
        elseif strcmp(behav_sequence{b_i}, 'R')
            behav_data = [behav_data, -1];
        end
    end
    cfg_switch = [];
    cfg_switch.devt_max = devt_max_behav;
    [p_switch_behav{sess_i}] = calculate_p_switch_by_index(cfg_switch, behav_data);
    % Calculate baseline: 2*pL*pR
    pL = sum(behav_data > 0) / length(behav_data);
    p_switch_baseline_behav(sess_i) = 2 * pL * (1-pL);
    
    adj_p_switch_behav(sess_i, :) = p_switch_behav{sess_i} - p_switch_baseline_behav(sess_i);
end

%% Simulate indepedent events for testing
n_evts = 1000;
pLs = [0.1, 0.2, 0.3, 0.4, 0.5];
out = cell(1, length(pLs));

for pL_i = 1:length(pLs)
    pL = pLs(pL_i);
    out{pL_i}.actual_pL = zeros(1, n_evts);
    out{pL_i}.actual_pR = zeros(1, n_evts);

    for evt_i = 1:n_evts
        if rand() <= pL
            out{pL_i}.actual_pL(evt_i) = 1;
        else
            out{pL_i}.actual_pR(evt_i) = 1;
        end
    end
end

%% Plotting probablity of switching with baseline for each session
for sess_i = 1:length(out)
    figure;
    plot(1:length(p_switch{sess_i}), p_switch{sess_i}, '.-r'); hold on;
    % Behavior
    plot(1:length(p_switch_behav{sess_i}), p_switch_behav{sess_i}, '.-g');

    u_bound = prctile(p_switch_shuffles{sess_i}, 97.5, 1) - mean(p_switch_shuffles{sess_i}, 1);
    l_bound = mean(p_switch_shuffles{sess_i}, 1) - prctile(p_switch_shuffles{sess_i}, 2.5, 1);
    h = shadedErrorBar(1:length(p_switch{sess_i}), mean(p_switch_shuffles{sess_i}, 1), [u_bound;l_bound]);
    yline(p_switch_baseline(sess_i), '--k', 'Baseline');

    ylim([0, 0.8]);
    xlabel('Delta event index')
    ylabel('P(interleave)')
    set(gca,'FontSize', 18)
    if sig_events_only
        title(sprintf('Significant events by odd ratio Session %d', sess_i));
    else
        title(sprintf('All events Session %d', sess_i));
    end
    % saveas(gcf, sprintf('sig_odds_%d.jpg', sess_i));
end

%%
adj_p_switch_shuffles_mat = zeros(1000, devt_max);
for i = 1:1000
    shuffles_mat = [];
    for sess_i = 1:length(out)
        if isempty(shuffles_mat)
            shuffles_mat = adj_p_switch_shuffles{sess_i}(i, :);
        else
            shuffles_mat = [shuffles_mat; adj_p_switch_shuffles{sess_i}(i, :)];
        end
    end
    adj_p_switch_shuffles_mat(i, :) = mean(shuffles_mat, 1);
end

% adj_p_switch_shuffles_m = cellfun(@(x) mean(x, 1) ,adj_p_switch_shuffles, 'UniformOutput', false);
%
% adj_p_switch_shuffles_mat = [];
% for sess_i = 1:length(out)
%     if isempty(adj_p_switch_shuffles_mat)
%         adj_p_switch_shuffles_mat = adj_p_switch_shuffles_m{sess_i};
%     else
%         adj_p_switch_shuffles_mat = [adj_p_switch_shuffles_mat; adj_p_switch_shuffles_m{sess_i}];
%     end
% end

%%
devt_max = 50;
devt_nbin = 5;

n_bins = devt_max / devt_nbin;
adj_p_switch_bin = zeros(size(adj_p_switch, 1), n_bins);
adj_p_switch_shuffles_mat_bin = zeros(size(adj_p_switch_shuffles_mat, 1), n_bins);

for bin_i = 0:n_bins-1
    devt_i_start = bin_i * devt_nbin + 1;
    devt_i_end = devt_i_start + devt_nbin - 1;
    adj_p_switch_bin(:, bin_i+1) = nanmean(adj_p_switch(:, devt_i_start:devt_i_end), 2);
    adj_p_switch_shuffles_mat_bin(:, bin_i+1) = nanmean(adj_p_switch_shuffles_mat(:, devt_i_start:devt_i_end), 2);
end

%% Plotting adjusted (by independent baseline) probablity of switching across sessions
adj_p_switch_m = nanmean(adj_p_switch_bin, 1);
adj_p_switch_sem = nanstd(adj_p_switch_bin, 1) / sqrt(size(adj_p_switch_bin, 1));

x = 1:size(adj_p_switch_bin, 2);
xpad = 0.5;
h = errorbar(x, adj_p_switch_m, adj_p_switch_sem, 'LineWidth', 1.5); hold on;
set(h, 'Color', 'k');

adj_p_switch_shuffles_m = mean(adj_p_switch_shuffles_mat_bin, 1);
u_bound = prctile(adj_p_switch_shuffles_mat_bin, 97.5, 1) - adj_p_switch_shuffles_m;
l_bound = adj_p_switch_shuffles_m - prctile(adj_p_switch_shuffles_mat_bin, 2.5, 1);
sh = shadedErrorBar(x, adj_p_switch_shuffles_m, [u_bound;l_bound]);

% sh = errorbar(x, shuf_adj_p_switch_m, shuf_adj_p_switch_sem, 'LineWidth', 1.5); hold on;
% set(sh, 'Color', [0.7 0.7 0.7]);

hold on;
plot(x, adj_p_switch_m, '.k', 'MarkerSize', 20);
set(gca, 'XTick', x, 'XTickLabel', 1:5:size(adj_p_switch, 2),  ...
    'XLim', [x(1)-xpad x(end)+xpad], 'YLim', [-0.03, 0.03], 'FontSize', 18, ...
    'LineWidth', 1, 'TickDir', 'out');
box off;
plot([x(1)-xpad devt_max+xpad], [0 0], '--k', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);

xlabel('Delta event index')
ylabel('adjusted P(interleave)')

if sig_events_only
    title('Significant events by odd ratio');
else
    title('All events');
end
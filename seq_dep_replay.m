%%
restoredefaultpath; % start with clean slate
addpath(genpath('/Users/mac/Projects/vandermeerlab/code-matlab/shared'));
addpath(genpath('/Users/mac/Projects/vandermeerlab/code-matlab/tasks/Alyssa_Tmaze'));
addpath(genpath('/Users/mac/Projects/papers-master/Carey_etal_submitted'));
addpath(genpath('/Users/mac/Projects/seq_dep_replay/utils'));
addpath(genpath('/Users/mac/Projects/seq_dep_replay/procedures'));

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
L_sig = out.shuf_perc >= 0.95 & out.shuf_perc < 1;
R_sig = out.shuf_perc <= 0.05 & out.shuf_perc > 0;

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

%% Plot Inter-SWRs-interval as a function of time.
SWR_data = actual_L_R_diff;
SWR_times = Q_SWR.tvec;

sig_SWR_idx = find(L_sig | R_sig);
sig_SWR_times = SWR_times(sig_SWR_idx);

cfg_sw = [];
% bin_size = 0.2;
% cfg_sw.bin_egdes = -1.2:bin_size:1.8;
% [p_switch, SWR_t_diffs] = calculate_p_switch(cfg_sw, SWR_data, SWR_times);

bin_size = 0.25;
cfg_sw.bin_egdes = 0:bin_size:2;
[p_switch, SWR_t_diffs] = calculate_p_switch(cfg_sw, SWR_data(sig_SWR_idx), sig_SWR_times);

t_diffs_x = cfg_sw.bin_egdes(1:end-1) + bin_size / 2;
plot(t_diffs_x, p_switch, '.-r');
yline(0.5, '--k', 'HandleVisibility','off');

xlim([-0.75, 1.75]);
ylim([-0.25, 1.25]);
xlabel('log10 (time elapsed since last SWR)')
ylabel('P(switch)')
set(gca,'FontSize', 18)

%% Shuffling SWR indices
n_shuffles = 1000;
p_switch_shuffles = zeros(n_shuffles, length(t_diffs_x));

for s_i = 1:n_shuffles
%     shuffle_indices = randperm(length(SWR_data));
%     [s_p_switch] = calculate_p_switch(SWR_data(shuffle_indices), SWR_times);

    shuffle_indices = randperm(length(sig_SWR_idx));
    s_sig_SWR_idx = sig_SWR_idx(shuffle_indices);
    % SWR times for significant events should stay constant
    [s_p_switch] = calculate_p_switch(cfg_sw, SWR_data(s_sig_SWR_idx), SWR_times(sig_SWR_idx));
    
    p_switch_shuffles(s_i, :) = s_p_switch;
end

u_bound = prctile(p_switch_shuffles, 97.5, 1) - mean(p_switch_shuffles, 1);
l_bound = mean(p_switch_shuffles, 1) - prctile(p_switch_shuffles, 2.5, 1);
h = shadedErrorBar(t_diffs_x, mean(p_switch_shuffles, 1), [u_bound;l_bound]);

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

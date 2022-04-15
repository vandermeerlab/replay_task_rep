%% Position (with spikes) plotting
figure;
for i = 1:length(pos{3})
    X = pos{3}{i}.data(:, 2);
    Y = pos{3}{i}.data(:, 3);
    
    keep_idx = X ~= 0 & Y ~= 0;
    subplot(2, 4, i);
    plot(X(keep_idx), Y(keep_idx), '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
    axis off; hold on;
    
    spike_times = spikes{3}{i}{end}{1}.data(:, 1);
    pos_times = pos{3}{i}.data(:, 1);
    spike_X = interp1(pos_times, X, spike_times, 'linear');
    spike_Y = interp1(pos_times, Y, spike_times, 'linear');
    
    plot(spike_X, spike_Y, '.r')
end

%% Extract trial information
session_i = 4;
pos_times = pos{3}{session_i}.data(:, 1);
X = pos{3}{session_i}.data(:, 2);
Y = pos{3}{session_i}.data(:, 3);

ExpKeys_filename = 'ExpKeys03.mat';
if isfile('ExpKeys03.mat')
    load(ExpKeys_filename);
else
    ExpKeys = cell(1, 7);
end

%% Get coordinates of center, left and right arm reward zones
plot(X, Y, '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;

[center_X, center_Y] = ginput(4);
plot(center_X, center_Y); hold on;

[left_X, left_Y] = ginput(4);
plot(left_X, left_Y); hold on;

[right_X, right_Y] = ginput(4);
plot(right_X, right_Y); hold on;

%%
ExpKeys{session_i}.center.x = [60, 90];
ExpKeys{session_i}.center.y = [123, 169];
ExpKeys{session_i}.left.x = [28, 54];
ExpKeys{session_i}.left.y = [123, 169];
ExpKeys{session_i}.right.x = [93, 130];
ExpKeys{session_i}.right.y = [123, 169];

%% Assign each position a label if in center: 1, left: 2, right: 4, otherwsie 0.
% Then we can use the difference between location labels to potential trial
% starts or ends. For example, -1 would be exiting center reward zone, etc.
pos_trial = zeros(1, length(X));
for pos_i = 1:length(pos_trial)
    % Check if position is in the center arm reward area
    if (X(pos_i) > ExpKeys{session_i}.center.x(1) && X(pos_i) < ExpKeys{session_i}.center.x(2))...
            && (Y(pos_i) > ExpKeys{session_i}.center.y(1) && Y(pos_i) < ExpKeys{session_i}.center.y(2))
        pos_trial(pos_i) = 1;
    % Check if position is in the left arm reward area
    elseif (X(pos_i) > ExpKeys{session_i}.left.x(1) && X(pos_i) < ExpKeys{session_i}.left.x(2))...
            && (Y(pos_i) > ExpKeys{session_i}.left.y(1) && Y(pos_i) < ExpKeys{session_i}.left.y(2))
        pos_trial(pos_i) = 2;
    % Check if position is in the right arm reward area
    elseif (X(pos_i) > ExpKeys{session_i}.right.x(1) && X(pos_i) < ExpKeys{session_i}.right.x(2))...
            && (Y(pos_i) > ExpKeys{session_i}.right.y(1) && Y(pos_i) < ExpKeys{session_i}.right.y(2))
        pos_trial(pos_i) = 4;
    else
        pos_trial(pos_i) = 0;
    end
end

diff_pos_trial = diff(pos_trial);

%% Locatons in the reward zone
plot(X, Y, '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;
% center
scatter(X(pos_trial == 1), Y(pos_trial == 1)); hold on;
% left
scatter(X(pos_trial == 2), Y(pos_trial == 2)); hold on;
% right
scatter(X(pos_trial == 4), Y(pos_trial == 4)); hold on;

%% Position of potential trial events
plot(X, Y, '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;
% center-exiting
scatter(X(diff_pos_trial == -1), Y(diff_pos_trial == -1)); hold on;
% left-entering
scatter(X(diff_pos_trial == 2), Y(diff_pos_trial == 2)); hold on;
% right-entering
scatter(X(diff_pos_trial == 4), Y(diff_pos_trial == 4)); hold on;

%% Extract left and right (outbound) running trials
cfg_outbound.mode = 'next';
this_session.center_exit = pos_times(diff_pos_trial == -1);
this_session.left_enter = pos_times(diff_pos_trial == 2);
this_session.right_enter = pos_times(diff_pos_trial == 4);

trial_start_L = [];
trial_end_L = [];
trial_start_R = [];
trial_end_R = [];

for c_i = 1:length(this_session.center_exit)
    [pos_t,fieldname] = FindFieldTime(cfg_outbound, this_session, this_session.center_exit(c_i));
    if ~strcmp(fieldname{1}, 'center_exit')
        if strcmp(fieldname{1}, 'left_enter')
            trial_start_L = [trial_start_L, this_session.center_exit(c_i)];
            trial_end_L = [trial_end_L, pos_t];
        elseif strcmp(fieldname{1}, 'right_enter')
            trial_start_R = [trial_start_R, this_session.center_exit(c_i)];
            trial_end_R = [trial_end_R, pos_t];
        end
    end
end

%%
ExpKeys{session_i}.trial_start_L = trial_start_L;
ExpKeys{session_i}.trial_end_L = trial_end_L;
ExpKeys{session_i}.trial_start_R = trial_start_R;
ExpKeys{session_i}.trial_end_R = trial_end_R;

%% Position of outbounding trials
figure;
all_trial_starts = sort([trial_start_L, trial_start_R]);
all_trial_ends = sort([trial_end_L, trial_end_R]);
for l_i = 1:length(all_trial_starts)
    keep_idx = (pos_times > all_trial_starts(l_i)) & (pos_times < all_trial_ends(l_i));
    subplot(3, 4, l_i);
    scatter(X(pos_times == all_trial_starts(l_i)), Y(pos_times == all_trial_starts(l_i)), 'g'); hold on;
    plot(X(keep_idx), Y(keep_idx), '.','Color',[0.5 0.5 0.5],'MarkerSize',1); hold on;
    scatter(X(pos_times == all_trial_ends(l_i)), Y(pos_times == all_trial_ends(l_i)), 'r'); hold on;
    axis off; hold on;
end

%%
save('ExpKeys03.mat','ExpKeys');
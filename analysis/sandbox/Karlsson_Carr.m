%% Position (with spikes) plotting
figure;
day_i = 2;
for i = 1:length(pos{day_i})
    X = pos{day_i}{i}.data(:, 2);
    Y = pos{day_i}{i}.data(:, 3);
    
    keep_idx = X ~= 0 & Y ~= 0;
    subplot(2, 4, i);
    plot(X(keep_idx), Y(keep_idx), '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
    axis off; hold on;
    
    spike_times = spikes{day_i}{i}{end}{1}.data(:, 1);
    pos_times = pos{day_i}{i}.data(:, 1);
    spike_X = interp1(pos_times, X, spike_times, 'linear');
    spike_Y = interp1(pos_times, Y, spike_times, 'linear');
    
    plot(spike_X, spike_Y, '.r')
end

%% Extract trial information
day_i = 2;
session_i = 4;
pos_times = pos{day_i}{session_i}.data(:, 1);
X = pos{day_i}{session_i}.data(:, 2);
Y = pos{day_i}{session_i}.data(:, 3);

[~, subj, ~] = fileparts(pwd);
ExpKeys_filename = sprintf('%sExpKeys0%d.mat', subj, day_i);
if isfile(ExpKeys_filename)
    load(ExpKeys_filename);
else
    ExpKeys = cell(1, length(pos{day_i}));
end

%% Get coordinates of center, left and right arm reward zones
figure;
plot(X, Y, '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;

% Label it in this order: top left -> bottom left -> bottom right -> top right
[center_X, center_Y] = ginput(4);
plot(center_X, center_Y); hold on;

[right_X, right_Y] = ginput(4);
plot(right_X, right_Y); hold on;

[left_X, left_Y] = ginput(4);
plot(left_X, left_Y); hold on;

%%
ExpKeys{session_i}.center.x = [floor(min(center_X(1:2))), ceil(max(center_X(3:4)))];
ExpKeys{session_i}.center.y = [floor(min(center_Y(2:3))), ceil(max(center_Y([1, 4])))];
ExpKeys{session_i}.right.x = [floor(min(right_X(1:2))), ceil(max(right_X(3:4)))];
ExpKeys{session_i}.right.y = [floor(min(right_Y(2:3))), ceil(max(right_Y([1, 4])))];
ExpKeys{session_i}.left.x = [floor(min(left_X(1:2))), ceil(max(left_X(3:4)))];
ExpKeys{session_i}.left.y = [floor(min(left_Y(2:3))), ceil(max(left_Y([1, 4])))];

%% Assign each position a label if in center: 1, left: 2, right: 4, otherwsie 0.
% Then we can use the difference between location labels to potential trial
% starts or ends. For example, -1 would be exiting center reward zone, etc.
pos_trial = zeros(1, length(X));
for pos_i = 1:length(pos_trial)
    % Check if position is in the center arm reward area
    if (X(pos_i) > ExpKeys{session_i}.center.x(1) && X(pos_i) < ExpKeys{session_i}.center.x(2))...
            && (Y(pos_i) > ExpKeys{session_i}.center.y(1) && Y(pos_i) < ExpKeys{session_i}.center.y(2))
        pos_trial(pos_i) = 1;
    % Check if position is in the right arm reward area
    elseif (X(pos_i) > ExpKeys{session_i}.right.x(1) && X(pos_i) < ExpKeys{session_i}.right.x(2))...
            && (Y(pos_i) > ExpKeys{session_i}.right.y(1) && Y(pos_i) < ExpKeys{session_i}.right.y(2))
        pos_trial(pos_i) = 4;
    % Check if position is in the left arm reward area
    elseif (X(pos_i) > ExpKeys{session_i}.left.x(1) && X(pos_i) < ExpKeys{session_i}.left.x(2))...
            && (Y(pos_i) > ExpKeys{session_i}.left.y(1) && Y(pos_i) < ExpKeys{session_i}.left.y(2))
        pos_trial(pos_i) = 2;
    else
        pos_trial(pos_i) = 0;
    end
end

diff_pos_trial = diff(pos_trial);

%% Locatons in the reward zone
figure;
plot(X, Y, '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;
% center
scatter(X(pos_trial == 1), Y(pos_trial == 1)); hold on;
% right
scatter(X(pos_trial == 4), Y(pos_trial == 4)); hold on;
% left
scatter(X(pos_trial == 2), Y(pos_trial == 2)); hold on;

%% Position of potential trial events
figure;
plot(X, Y, '.','Color',[0.5 0.5 0.5],'MarkerSize',1);
axis off; hold on;
% center-exiting
scatter(X(diff_pos_trial == -1), Y(diff_pos_trial == -1)); hold on;
% right-entering
scatter(X(diff_pos_trial == 4), Y(diff_pos_trial == 4)); hold on;
% left-entering
scatter(X(diff_pos_trial == 2), Y(diff_pos_trial == 2)); hold on;

%% Extract left and right (outbound) running trials
cfg_outbound.mode = 'next';
this_session = {};
this_session.center_exit = pos_times(diff_pos_trial == -1);
this_session.right_enter = pos_times(diff_pos_trial == 4);
this_session.left_enter = pos_times(diff_pos_trial == 2);

trial_start_R = [];
trial_end_R = [];
trial_start_L = [];
trial_end_L = [];

for c_i = 1:length(this_session.center_exit)
    [pos_t,fieldname] = FindFieldTime(cfg_outbound, this_session, this_session.center_exit(c_i));
    if ~isempty(fieldname) && ~strcmp(fieldname{1}, 'center_exit')
        if strcmp(fieldname{1}, 'right_enter')
            trial_start_R = [trial_start_R, this_session.center_exit(c_i)];
            trial_end_R = [trial_end_R, pos_t];
        elseif strcmp(fieldname{1}, 'left_enter')
            trial_start_L = [trial_start_L, this_session.center_exit(c_i)];
            trial_end_L = [trial_end_L, pos_t];
        end
    end
end

%% Position of outbounding trials
figure;
all_trial_starts = sort([trial_start_R, trial_start_L]);
all_trial_ends = sort([trial_end_R, trial_end_L]);
for l_i = 1:length(all_trial_starts)
    keep_idx = (pos_times > all_trial_starts(l_i)) & (pos_times < all_trial_ends(l_i));
    subplot(ceil(length(all_trial_starts) / 4), 4, l_i);
    scatter(X(pos_times == all_trial_starts(l_i)), Y(pos_times == all_trial_starts(l_i)), 'g'); hold on;
    plot(X(keep_idx), Y(keep_idx), '.','Color',[0.5 0.5 0.5],'MarkerSize',1); hold on;
    scatter(X(pos_times == all_trial_ends(l_i)), Y(pos_times == all_trial_ends(l_i)), 'r'); hold on;
    axis off; hold on;
    title(sprintf('duration: %.2f s', all_trial_ends(l_i) - all_trial_starts(l_i)));
end

%%
ExpKeys{session_i}.trial_start_R = trial_start_R;
ExpKeys{session_i}.trial_end_R = trial_end_R;
ExpKeys{session_i}.trial_start_L = trial_start_L;
ExpKeys{session_i}.trial_end_L = trial_end_L;

%%
save(ExpKeys_filename,'ExpKeys');


%% Count the number of CA1 cells for all days across subjects.
temp_fd = dir;
temp_fd = temp_fd(3:end);
temp_fd = temp_fd([temp_fd.isdir]);
n_CA1_list = [];

for iFD = 1:length(temp_fd)
    cd(temp_fd(iFD).name)
    m = FindFiles('*cellinfo.mat');
    if ~isempty(m)
        load(m{1})
    end
    for day_i = 1:length(cellinfo)
        if ~isempty(cellinfo{day_i})
            all_tetrodes = cellinfo{day_i}{2};
            n_CA1_cells = 0;
            for tet_i = 1:length(all_tetrodes)
                if ~isempty(all_tetrodes{tet_i})
                    this_tetrode = all_tetrodes{tet_i};
                    for cell_i = 1:length(this_tetrode)
                        if ~isempty(this_tetrode{cell_i}) && isfield(this_tetrode{cell_i}, 'area') && strcmp(this_tetrode{cell_i}.area, 'CA1')
                            n_CA1_cells = n_CA1_cells + 1;
                        end
                    end
                end
            end
            if n_CA1_cells >= 30
                disp(temp_fd(iFD).name)
                disp(day_i)
            end
            n_CA1_list = [n_CA1_list, n_CA1_cells];
        end
    end
    clear m
    cd ..
end

%% Count the number of sessions whose number of CA1 cells is larger than threshold (40).

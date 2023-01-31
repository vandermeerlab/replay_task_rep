LoadExpKeys;

evt = LoadEvents([]);
% Quick check to remove empty label
non_empty_idx = ~cellfun(@isempty, evt.label);
evt.label = evt.label(non_empty_idx);

please = []; please.load_questionable_cells = 1; please.getTTnumbers = 1;
S = LoadSpikes(please);

if isfield(ExpKeys,'TetrodeTargets')
    % keep only hippocampus cells
    hc_tt = find(strcmp(ExpKeys.Target, 'Hippocampus'));
    hc_tt = find(ExpKeys.TetrodeTargets == hc_tt);
    keep_idx = ismember(S.usr.tt_num, hc_tt);
    S = SelectTS([], S, keep_idx);
end

%% find left/right (rewarded) trial times
%keep = ~cellfun('isempty',evt.label); evt = SelectTS([],evt,keep); % may fail with current codebase master, works with striatal-spike-rhythms repo

min_trial_len = 1; % in seconds, used to remove multiple feeder fires
if isfield(ExpKeys,'FeederL1') % feeder IDs defined, use them

    feeders = cat(2, ExpKeys.FeederL1, ExpKeys.FeederR1);
    feeder_labels = {'L', 'R'};
    reward_t = [];
    ll = @(x) x(end); % function to get last character of input
    for iF = 1:length(feeders)

        keep_idx = find(num2str(feeders(iF)) == cellfun(ll, evt.label));
        reward_t.(feeder_labels{iF}) = evt.t{keep_idx};

        % remove multiple feeder fires
        ifi = cat(2, Inf, diff(reward_t.(feeder_labels{iF})));
        reward_t.(feeder_labels{iF}) = reward_t.(feeder_labels{iF})(ifi >= min_trial_len);

    end

else
    error('no left/right feeder IDs defined');
end

%% find the trial order to look like MotivationalT metadata
nL = length(reward_t.L); nR = length(reward_t.R);
left_labels = repmat({'L'}, [1 nL]); right_labels = repmat({'R'}, [1 nR]);
all_labels = cat(2, left_labels, right_labels);
all_times = cat(2, reward_t.L, reward_t.R);

[sorted_reward_times, sort_idx] = sort(all_times, 'ascend');
sequence = all_labels(sort_idx);

%% convert reward times into trial ivs
trial_len = 2; % start and end are these many seconds relaive to feeder fire

trial_iv_L = iv(reward_t.L - trial_len, reward_t.L);
trial_iv_R = iv(reward_t.R - trial_len, reward_t.R);

trial_iv = iv(all_times - trial_len, all_times);

%% should now be able to use GetMatchedTrials()
metadata = [];
metadata.taskvars.trial_iv = trial_iv;
metadata.taskvars.trial_iv_L = trial_iv_L;
metadata.taskvars.trial_iv_R = trial_iv_R;
metadata.taskvars.sequence = sequence;

ExpKeys.badTrials = [];

%% Plot some example trajectories
pos = LoadPos([]);
pos_X = getd(pos,'x');
pos_Y = getd(pos,'y');

figure;
plot(pos_X, pos_Y);
hold on;

%% Get MS (trial start) position
[MS_X, MS_Y] = ginput(4);
plot(MS_X, MS_Y); hold on;

%%
ExpKeys.MS.x = [floor(min(MS_X(1:2))), ceil(max(MS_X(3:4)))];
ExpKeys.MS.y = [floor(min(MS_Y(2:3))), ceil(max(MS_Y([1, 4])))];

save(FindFile('*keys.m'),'ExpKeys');

%% Assign each position a label as 1 if in MS
% Then we can use the difference between location labels to potential trial starts.
% For example, -1 would be exiting MS.
pos_trial = zeros(1, length(pos.data));
for pos_i = 1:length(pos_trial)
    if (pos_X(pos_i) > ExpKeys.MS.x(1) && pos_X(pos_i) < ExpKeys.MS.x(2))...
        && (pos_Y(pos_i) > ExpKeys.MS.y(1) && pos_Y(pos_i) < ExpKeys.MS.y(2))
        pos_trial(pos_i) = 1;
    end
end
diff_pos_trial = diff(pos_trial);

%% Locatons in the MS zone
scatter(pos_X(pos_trial == 1), pos_Y(pos_trial == 1)); hold on;

%% Position of potential start events
scatter(pos_X(diff_pos_trial == -1), pos_Y(diff_pos_trial == -1)); hold on;

%% Extract left and right running trial starts
cfg_trial.mode = 'prev';
this_session = {};
this_session.reward_times = all_times;
this_session.MS_exit = pos.tvec(diff_pos_trial == -1);

trial_starts = [];

for c_i = 1:length(this_session.reward_times)
    [pos_t,fieldname] = FindFieldTime(cfg_trial, this_session, this_session.reward_times(c_i));
    if ~isempty(fieldname) && ~strcmp(fieldname{1}, 'reward_times')
        if strcmp(fieldname{1}, 'MS_exit')
            trial_starts = [trial_starts, pos_t];
        end
    end
end

L_idx = cellfun(@(x) strcmp(x, 'L'), all_labels);
R_idx = cellfun(@(x) strcmp(x, 'R'), all_labels);

%% 
res_pos = restrict(pos,metadata.taskvars.trial_iv); % restricted interval only
plot(getd(res_pos,'x'),getd(res_pos,'y'),'r.');
% 
% sess_path = strsplit(pwd, '\');
% sess_name = sess_path{end};
% saveas(gcf, fullfile('C:\Users\mvdmlab\Desktop\adr_trajectories', sess_name), 'jpeg');

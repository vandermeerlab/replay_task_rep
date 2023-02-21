LoadExpKeys;
evt = LoadEvents([]);

load(FindFile('*Metadata.mat'))
metadata = Metadata;
clear Metadata;

%% Quick check to remove empty label
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
    % Some oddity in Guptat dataset, Feeder labels are inversed.
    feeder_labels = {'R', 'L'};
%     feeder_labels = {'L', 'R'};
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

%% Plot some example trajectories
pos = LoadPos([]);
pos_X = getd(pos,'x');
pos_Y = getd(pos,'y');

figure;
plot(pos_X, pos_Y,'.','Color',[0.7 0.7 0.7],'MarkerSize',4); xlabel('x data'); ylabel('y data');
view(90,90);
hold on;

%% Get MS (trial start) position
[MS_X, MS_Y] = ginput(4);
plot(MS_X, MS_Y); hold on;

%%
Metadata.MS.x = [floor(min(MS_X(1:2))), ceil(max(MS_X(3:4)))];
Metadata.MS.y = [floor(min(MS_Y(2:3))), ceil(max(MS_Y([1, 4])))];

%% Assign each position a label as 1 if in MS
% Then we can use the difference between location labels to potential trial starts.
% For example, -1 would be exiting MS.
pos_trial = zeros(1, length(pos.data));
for pos_i = 1:length(pos_trial)
    if (pos_X(pos_i) > Metadata.MS.x(1) && pos_X(pos_i) < Metadata.MS.x(2))...
        && (pos_Y(pos_i) > Metadata.MS.y(1) && pos_Y(pos_i) < Metadata.MS.y(2))
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
this_session.reward_times = sorted_reward_times;
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

L_idx = cellfun(@(x) strcmp(x, 'L'), sequence);
R_idx = cellfun(@(x) strcmp(x, 'R'), sequence);

%% Plot trial data
trial_iv = iv(trial_starts, sorted_reward_times);
trial_iv_L = iv(trial_starts(L_idx), sorted_reward_times(L_idx));
trial_iv_R = iv(trial_starts(R_idx), sorted_reward_times(R_idx));

res_pos = restrict(pos, trial_iv_L); % restricted interval only
plot(getd(res_pos,'x'),getd(res_pos,'y'),'r.');

%% should now be able to use GetMatchedTrials()
Metadata.taskvars.trial_iv = trial_iv;
Metadata.taskvars.trial_iv_L = trial_iv_L;
Metadata.taskvars.trial_iv_R = trial_iv_R;
Metadata.taskvars.sequence = sequence;

%% Now get the conversion factors using the function PosCon()

% PosCon() will reorient your maze if you pass in certain varargins. 
% read the help documentation on PosCon for more intructions

realTrackDims = [108 160]; % x width and y width in centimeters
convFact = PosCon(pos,realTrackDims, 'rot', 90);

% You should save convFact in ExpKeys
%          ExpKeys.convFact = [xConvFact yConvFact];

%% Load position data in units of centimeters
% Now that we have the conversion factor for pixels to centimeters, we can also plot data
% in cm using the .convFact flag in LoadPos()
cfg = [];
cfg.convFact = convFact;

pos_cm = LoadPos(cfg);
figure; plot(getd(pos_cm,'x'),getd(pos_cm,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4);
xlabel('x data'); ylabel('y data');
view(90,90);
title('Position in centimeters');
hold on;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                    Making a Coord File                              %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open the function MakeCoord() and read its internal documentation for
% more information

% let's get the idealized trajectories; MakeCoord() takes varargins that can
% reorient your maze as you prefer to see it.
coordL = MakeCoord(pos,'titl','Draw left trajectory, press enter when done', 'rot', 90); % CoordL is in units of pixels
coordR = MakeCoord(pos,'titl','Draw right trajectory, press enter when done', 'rot', 90); % CoordR is in units of pixels


%% plot the coords
figure('units','normalized','outerposition',[0,0.25,1,.75]); %a nice way of specifying figure size and position
subplot(121)
hold on;
plot(pos,'.','MarkerSize',4,'Color',[.7 .7 .7]);
plot(coordL.coord(1,:),coordL.coord(2,:),'o','MarkerSize',3);
xlabel('x data'); ylabel('y data'); title('Left coords in px');
view(90,90);

subplot(122)
hold on;
plot(pos,'.','MarkerSize',4,'Color',[.7 .7 .7]);
plot(coordR.coord(1,:),coordR.coord(2,:),'o','MarkerSize',3);
xlabel('x data'); ylabel('y data'); title('Right coords in px');
view(90,90);

%% conver coords into cm
% these coords should be converted to units of cm using the convFact you already
% collected above:

coordL_cm = coordL; % copy coordL under a new variable name, and apply some changes:
coordL_cm.coord(1,:) = coordL_cm.coord(1,:)./convFact(1); % apply x conversion
coordL_cm.coord(2,:) = coordL_cm.coord(2,:)./convFact(2); % apply y conversion
coordL_cm.units = 'cm';

coordR_cm = coordR; % as above, for R instead
coordR_cm.coord(1,:) = coordR_cm.coord(1,:)./convFact(1); % apply x conversion
coordR_cm.coord(2,:) = coordR_cm.coord(2,:)./convFact(2); % apply y conversion
coordR_cm.units = 'cm';

% put it all in a struct for tighter packing in the base workspace (when loading variables later)
coord = struct('coordL',coordL,'coordL_cm',coordL_cm,'coordR',coordR,'coordR_cm',coordR_cm);

clear coordL coordL_cm coordR coordR_cm

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                   Getting Choice Points Manually                    %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If your maze has any choice points, it's probably a good idea to get the
% coordinates of those choice points. Here's how you can do this using a
% script:

figure; plot(getd(pos,'x'),getd(pos,'y'),'.','Color',[0.7 0.7 0.7],'MarkerSize',4); 
view(90,90);
hold on;
plot(coord.coordL.coord(1,:),coord.coordL.coord(2,:),'ob'); 
plot(coord.coordR.coord(1,:),coord.coordR.coord(2,:),'og'); title('Click choice point; press enter');
maximize

% get user input:
[x,y] = ginput;

plot(x,y,'or','MarkerSize',10,'LineWidth',4); pause(1); close

% convert choice point units

chp = [x; y];
chp_cm = [x/convFact(1); y/convFact(2)];

% add to coord
coord.chp = chp;
coord.chp_cm = chp_cm;

% coord can be saved as a metadata field, and metadata can be saved as a
% .mat file for later use

%%
save(FindFile('*metadata.mat'), 'metadata')

%% Plot some example trajectories
pos = LoadPos([]);
pos_X = getd(pos,'x');
pos_Y = getd(pos,'y');

figure;
plot(pos_X, pos_Y,'.','Color',[0.7 0.7 0.7],'MarkerSize',4); xlabel('x data'); ylabel('y data');
view(90,90);
hold on;

%% Plot linearized position with trial data
res_pos = restrict(pos, Metadata.taskvars.trial_iv_R); % restricted interval only
res_coord = coord.coordR;
% use LinearizePos() (NOTE: both our position tsd and coords are in cm!)
cfg = [];
cfg.debugMode = 1;
[linpos] = LinearizePos(cfg,res_pos,res_coord);

linpos.data = linpos.data(2, :);
n_all_samples = size(linpos.data, 2);

% Exclude position data that deviate from linearized path too much
cfg_path = []; cfg_path.method = 'raw'; cfg_path.threshold = 20; cfg_path.operation = '<';
path_iv = TSDtoIV(cfg_path, linpos);

[linpos] = LinearizePos([],res_pos,res_coord);

linpos = restrict(linpos, path_iv);
n_dev_samples = size(linpos.data, 2);
res_pos = restrict(res_pos, path_iv);

scatter(getd(res_pos,'x'),getd(res_pos,'y'), [], linpos.data); colorbar;
title(sprintf('%d out of %d samples (%.2f%%) remaining with %d pixels away included.', ...
    n_dev_samples, n_all_samples, n_dev_samples*100/n_all_samples, cfg_path.threshold))
hold on;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%                        Standardizing coords                         %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% StandardizeCoord takes a raw coord file and returns a standardized coord. 
% To do this we need the true path length of the used defined trajectory.
% This value is usually measured by the experimented and stored in the 
% ExpKeys.

run_dist = 250;

% use StandardizeCoord()
cfg = [];
coord_std = StandardizeCoord(cfg,coord.coordR_cm,run_dist);

% we can also specify arguments for how we want the coord to be standardized
coord_std2 = StandardizeCoord(cfg,coord.coordR_cm,run_dist,'pointDist',3);

%% make linearized position with standardized coords 
cfg = [];
linpos_std = LinearizePos(cfg,res_pos_cm,coord_std);
figure; plot(linpos_std,'.')

% compare with "raw" coords
lp1 = restrict(linpos,Metadata.taskvars.trial_iv_R.tstart(3),Metadata.taskvars.trial_iv_R.tend(3));
lp2 = restrict(linpos_std,Metadata.taskvars.trial_iv_R.tstart(3),Metadata.taskvars.trial_iv_R.tend(3));

figure;
hold on;
subplot(121); plot(lp1,'.'); axis tight; title('Linearized with raw coord');
% You can really tell that having less coord points "bins" the data. Depending on whether 
xlabel('Time (sec)'); ylabel('Position (idx)');
subplot(122); plot(lp2,'.'); axis tight; title('Linearized with standardized coord');
xlabel('Time (sec)'); ylabel('Position (idx)');
maximize
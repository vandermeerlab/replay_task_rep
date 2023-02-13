function [TC] = get_tuning_curve(cfg_in, session_path)
    % Adapted from https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/example_workflows/WORKFLOW_PlotOrderedRaster.m
    
    cfg_def.use_Gupta_data = 0;
    cfg_def.use_matched_trials = 1;
    cfg_def.removeInterneurons = 0;
    cfg_def.int_thres = 10;
    cfg_def.minSpikes = 25;

    mfun = mfilename;
    cfg = ProcessConfig(cfg_def,cfg_in,mfun);

    % Get the data
    cd(session_path);
    LoadExpKeys();
    if cfg.use_Gupta_data
        load(FindFile('*Metadata.mat'));
        metadata = Metadata;
        ExpKeys.badTrials = [];
    else
        LoadMetadata();
    end
    pos = LoadPos([]);

    cfg_spikes = {};
    cfg_spikes.load_questionable_cells = 1;
    S = LoadSpikes(cfg_spikes);

    if cfg.removeInterneurons
        channels = FindFiles('*.Ncs');
        cfg_lfp = []; cfg_lfp.fc = {channels{1}};
        lfp = LoadCSC(cfg_lfp);

        cfg_int = []; cfg_int.showFRhist = 0;
        cfg_int.max_fr = cfg.int_thres;
        S = RemoveInterneuronsHC(cfg_int,S, lfp);
    end

    %% set up data structs for 2 experimental conditions -- see lab wiki for this task at:
    % http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:task:motivationalt
    clear expCond;

    expCond(1).label = 'left'; % this is a T-maze, we are interested in 'left' and 'right' trials
    expCond(2).label = 'right'; % these are just labels we can make up here to keep track of which condition means what

    if cfg.use_matched_trials
        [matched_left, matched_right] = GetMatchedTrials({}, metadata, ExpKeys);
        expCond(1).t = matched_left;
        expCond(2).t = matched_right;
        tsxtart = [matched_left.tstart; matched_right.tstart];
        tend = [matched_left.tend; matched_right.tend];
    else
        expCond(1).t = metadata.taskvars.trial_iv_L; % previously stored trial start and end times for left trials
        expCond(2).t = metadata.taskvars.trial_iv_R;
        tstart = metadata.taskvars.trial_iv.tstart;
        tend = metadata.taskvars.trial_iv.tend;
    end

    expCond(1).coord = metadata.coord.coordL; % previously user input idealized linear track
    expCond(2).coord = metadata.coord.coordR; % note, this is in units of "camera pixels", not cm

    % % Remove cells with insufficient spikes
    % S_matched = restrict(S, tstart, tend);

    % spk_count = getSpikeCount([], S_matched);
    % cell_keep_idx = spk_count >= cfg.minSpikes;

    % S = SelectTS([], S, cell_keep_idx);

    expCond(1).S = S;
    expCond(2).S = S;

    %% linearize paths (snap x,y position samples to nearest point on experimenter-drawn idealized track)
    nCond = length(expCond);
    for iCond = 1:nCond
        this_coord.coord = expCond(iCond).coord; this_coord.units = 'px'; this_coord.standardized = 0;
        expCond(iCond).linpos = LinearizePos([],pos,this_coord.coord);
        % Compute the coordinate that the choice point corresponds
        chp = tsd(0,metadata.coord.chp,{'x','y'});
        chp.units = 'px';
        expCond(iCond).cp = LinearizePos([],chp,this_coord.coord);
    end

    %% Only include position data that not deviate from linearized path too much
    for iCond = 1:nCond
        this_coord.coord = expCond(iCond).coord; this_coord.units = 'px'; this_coord.standardized = 0;
        cfg_pos = [];
        cfg_pos.debugMode = 1;
        [linpos_path] = LinearizePos(cfg_pos,pos,this_coord.coord);
        % Position of deviation from linearized path
        linpos_path.data = linpos_path.data(2, :);
        cfg_path = []; cfg_path.method = 'raw'; cfg_path.threshold = 20; cfg_path.operation = '<';
        expCond(iCond).path_iv = TSDtoIV(cfg_path, linpos_path);
    end

    %% find intervals where rat is running
    spd = getLinSpd([],pos); % get speed (in "camera pixels per second")

    cfg_spd = []; cfg_spd.method = 'raw'; cfg_spd.threshold = 15;
    run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above 5 pix/s

    %% restrict (linearized) position data and spike data to desired intervals
    for iCond = 1:nCond
        
        fh = @(x) restrict(x,expCond(iCond).path_iv); % restrict S and linpos to reasonable path only
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});

        fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});

        fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    end

    %% get tuning curves, see lab wiki at:
    % http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week12
    for iCond = 1:nCond
        cfg_tc = []; cfg_tc.smoothingKernel = gausskernel(11, 1); cfg_tc.minOcc = 0.25;
        expCond(iCond).tc = TuningCurves(cfg_tc,expCond(iCond).S,expCond(iCond).linpos);
        % Temporal fix
        expCond(iCond).tc.tc(isnan(expCond(iCond).tc.tc)) = 0;
        [~,expCond(iCond).cp_bin] = histc(expCond(iCond).cp.data, expCond(iCond).tc.usr.binEdges);

    end
    
    for iCond = 1:nCond
        TC.(expCond(iCond).label).tc = expCond(iCond).tc.tc;
        TC.(expCond(iCond).label).cp_bin = expCond(iCond).cp_bin;
    end
end

% %% Plot animal's positions where a given neurons fires each spike.
% iC = 59;
% 
% cfg_spikes = {};
% cfg_spikes.load_questionable_cells = 1;
% S = LoadSpikes(cfg_spikes);
% S = restrict(S, Metadata.taskvars.trial_iv);
% S = restrict(S, run_iv);
% 
% spk_x = interp1(pos.tvec,getd(pos,'x'),S.t{iC},'linear');
% spk_y = interp1(pos.tvec,getd(pos,'y'),S.t{iC},'linear');
%  
% h = plot(spk_x,spk_y,'.r');

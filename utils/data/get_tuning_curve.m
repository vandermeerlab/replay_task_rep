function [TC] = get_tuning_curve(cfg_in, session_path)
    % Adapted from https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/example_workflows/WORKFLOW_PlotOrderedRaster.m
    
    cfg_def.use_Gupta_data = 0;
    cfg_def.use_matched_trials = 1;
    cfg_def.removeInterneurons = 0;
    cfg_def.minSpikes = 25;
    cfg_def.trackExcludeStart = 0; % exclude this amount (in px) from start and end of track
    cfg_def.trackExcludeEnd = 25;

    mfun = mfilename;
    cfg = ProcessConfig(cfg_def,cfg_in,mfun);

    % Get the data
    cd(session_path);
    LoadExpKeys();
    LoadMetadata();
    pos = LoadPos([]);

    cfg_spikes = {};
    cfg_spikes.load_questionable_cells = 1;
    S = LoadSpikes(cfg_spikes);

    if cfg.removeInterneurons
        cfg_lfp = []; cfg_lfp.fc = ExpKeys.goodSWR(1);
        lfp = LoadCSC(cfg_lfp);

%         cfg_int = []; cfg_int.showFRhist = 0;
%         cfg_int.max_fr = cfg.int_thres;
        S = RemoveInterneuronsHC([], S, lfp);
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
        tstart = [matched_left.tstart; matched_right.tstart];
        tend = [matched_left.tend; matched_right.tend];
    else
        expCond(1).t = metadata.taskvars.trial_iv_L; % previously stored trial start and end times for left trials
        expCond(2).t = metadata.taskvars.trial_iv_R;
        tstart = metadata.taskvars.trial_iv.tstart;
        tend = metadata.taskvars.trial_iv.tend;
    end

    expCond(1).coord = metadata.coord.coordL; % previously user input idealized linear track
    expCond(2).coord = metadata.coord.coordR; % note, this is in units of "camera pixels", not cm

    % remove cells with insufficient spikes
    [S] = removeEmptyCells(S);
    spk_count = getSpikeCount([], S);
    cell_keep_idx = spk_count >= cfg.minSpikes;
    S = SelectTS([], S, cell_keep_idx);

    expCond(1).S = S;
    expCond(2).S = S;
    expComb.S = S;

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
    run_iv = TSDtoIV(cfg_spd,spd); % intervals with speed above specified px/s

    %% exclude positions at beginning and end of track
    cfg_track1 = []; cfg_track1.method = 'raw'; cfg_track1.threshold = cfg.trackExcludeStart;
    cfg_track2 = []; cfg_track2.method = 'raw'; cfg_track2.operation = '<';

    for iCond = 1:nCond
        expCond(iCond).track_iv1 = TSDtoIV(cfg_track1,expCond(iCond).linpos);
        cfg_track2.threshold = max(expCond(iCond).linpos.data) - cfg.trackExcludeEnd;
        expCond(iCond).track_iv2 = TSDtoIV(cfg_track2,expCond(iCond).linpos);
    end

    %% restrict (linearized) position data and spike data to desired intervals
    for iCond = 1:nCond
        fh = @(x) restrict(x, expCond(iCond).track_iv1); % restrict S and linpos with beginning and end of track excluded
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});

        fh = @(x) restrict(x, expCond(iCond).track_iv2); % restrict S and linpos with beginning and end of track excluded
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});

        fh = @(x) restrict(x,expCond(iCond).path_iv); % restrict S and linpos to reasonable path only
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});

        fh = @(x) restrict(x,run_iv); % restrict S and linpos to run times only
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});

        fh = @(x) restrict(x,expCond(iCond).t); % restrict S and linpos to specific trials (left/right)
        expCond(iCond) = structfunS(fh,expCond(iCond),{'S','linpos'});
    end
    
    left_track_IV = IntersectIV([], expCond(1).track_iv2, expCond(1).track_iv1);
    right_track_IV = IntersectIV([], expCond(2).track_iv2, expCond(2).track_iv1);

    expComb.S = restrict(expComb.S, UnionIV([], expCond(1).path_iv, expCond(2).path_iv));
    expComb.S = restrict(expComb.S, run_iv);
    expComb.S = restrict(expComb.S, UnionIV([],left_track_IV, right_track_IV));
    expComb.t = UnionIV([],expCond(1).t,expCond(2).t);
    expComb.S = restrict(expComb.S, expComb.t);

    %% get tuning curves, see lab wiki at:
    % http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:nsb2015:week12
    for iCond = 1:nCond
        cfg_tc = []; cfg_tc.smoothingKernel = gausskernel(11, 1); cfg_tc.minOcc = 0.25;
        expCond(iCond).tc = TuningCurves(cfg_tc,expCond(iCond).S,expCond(iCond).linpos);
        % Temporal fix
        expCond(iCond).tc.tc(isnan(expCond(iCond).tc.tc)) = 0;
        [~,expCond(iCond).cp_bin] = histc(expCond(iCond).cp.data, expCond(iCond).tc.usr.binEdges);
        
        TC.(expCond(iCond).label) = expCond(iCond);
    end
    %%
    temp = expCond(2).linpos;
    temp.data = temp.data + max(expCond(1).linpos.data);
    % max(enc_TC.left.linpos.data);
    expComb.linpos = UnionTSD([], expCond(1).linpos, temp);
    expComb.linpos = restrict(expComb.linpos, expComb.t);

    mn = min(expComb.linpos.data);
    mx = max(expComb.linpos.data);
    expComb.nBins = 200;
    cfg_tc.binEdges{1} = linspace(mn,mx, expComb.nBins+1);

    TC.combined = expComb;
    TC.combined.tc = TuningCurves(cfg_tc, expComb.S, expComb.linpos);
end
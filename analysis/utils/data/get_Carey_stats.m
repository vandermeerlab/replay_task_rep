rng('default');

%%
fd = sort(getTmazeDataPath([]));

%% L vs R classification of run data (Figure S2b)
% first, grab relevant fields from output structure
cfg = []; cfg.sess = [2:6 10:11 13:24];
cfg.fn = fieldnames(out{1}); cfg.fn = cfg.fn(6:17); % relevant fields are the first 12 names

for iF = 1:length(cfg.fn)
   
    for iSess = 1:length(cfg.sess)
        this_sess = cfg.sess(iSess);
        
        eval(cat(2,cfg.fn{iF},'(',num2str(iSess),') = out{',num2str(this_sess),'}.',cfg.fn{iF},';'));
        
    end
    
end

xd = [1 2 4 5];
ywhat = {'left_preCP_correct','left_postCP_correct','right_preCP_correct','right_postCP_correct'};
ychance = {'left_chance_pre','left_chance_post','right_chance_pre','right_chance_post'};
cols = 'rrbb';

figure;
subplot(231);

for iW = 1:length(ywhat)
    
   this_ydata = eval(ywhat{iW}); this_ychance = eval(ychance{iW});
   h = bar(xd(iW),nanmean(this_ydata)); set(h,'FaceColor',cols(iW),'EdgeColor','none','BarWidth',0.6);
   hold on;
   h2 = plot(repmat(xd(iW),[1 length(this_ydata)]),this_ydata,'.k');
   h3 = plot([xd(iW)-0.4 xd(iW)+0.4],[nanmean(this_ychance) nanmean(this_ychance)],':k','LineWidth',2);
    
   % report means and SEMs; statistics vs. chance
   m = nanmean(this_ydata); s = nanstd(this_ydata)./sqrt(4); p = ranksum(this_ydata,this_ychance);
   fprintf('%s: %.2f +/- %.2f (SEM), p = %2.2e vs. chance\n',ywhat{iW},m,s,p);
   
end
set(gca,'LineWidth',1,'TickDir','out','FontSize',18,'XLim',[0 6],'XTick',xd,'XTickLabel',{'preL','postL','preR','postR'},'YLim',[0 1],'YTick',0:0.25:1);
box off;

%% collect data for L vs R SWR decoding
cfg = [];
cfg.sess = [2:6 10:11 13:24];
cfg.cutoff = 0.05; % percentile to include

% initialize variables to append to later
clear data;
rats = {'all','R042','R044','R050','R064'};
what = {'all','pre','task','post'};
vars = {'fracL_evt','fracR_evt','median_z','median_perc', ...
    'this_sess','this_trials','this_choice','this_type'}; % type should be last
for iR = 1:length(rats)
   for iW = 1:length(what)
       for iV = 1:length(vars)
    
       eval(sprintf('data.%s.%s.%s = [];',what{iW},rats{iR},vars{iV}));
       
       end
   end
end

% loop over sessions to populate variables
for iSess = 1:length(cfg.sess)
    
    this_sess = cfg.sess(iSess);
    cd(fd{this_sess}); LoadExpKeys; LoadMetadata;
    
    switch ExpKeys.RestrictionType
        case 'water'
            this_type = 1; % remember, water is on the right
        case 'food'
            this_type = 0; % food on the left
    end
    
    % behavior
    this_trials = length(strmatch('L',metadata.taskvars.sequence));
    this_trials = this_trials./length(metadata.taskvars.sequence);
  
    choice_trial_idx = setdiff(1:length(metadata.taskvars.sequence),ExpKeys.forcedTrials);
    choice_trials = metadata.taskvars.sequence(choice_trial_idx);
    this_choice = length(strmatch('L',choice_trials))./length(choice_trials);
    
    this_data = out{this_sess};
    this_rat = ExpKeys.goodSWR{1}(1:4);
    
    evt_times = ts;%({this_data.tvec});
    evt_times.t{1} = this_data.tvec;
    
    fprintf('Entering session %d...\n',this_sess);
    
    for iW = 1:length(what)
        
        fprintf('.%s\n',what{iW});
        
        % restrict to this epoch of interest
        clear keep;
        switch what{iW}
            case 'all'
                % do nothing, keep all data
                keep{1} = 1:length(this_data.tvec);
            case 'pre'
                [~,keep] = restrict_idx(evt_times,0,ExpKeys.prerecord(2));
            case 'task'
                [~,keep] = restrict_idx(evt_times,metadata.taskvars.rest_iv);
            case 'post'
                [~,keep] = restrict_idx(evt_times,ExpKeys.postrecord(1),Inf);
        end
        
        this_dataR = this_data;
        fnames = {'shuf_perc_odds','shuf_z_odds','tvec'};
        for iF = 1:length(fnames)
            temp = this_dataR.(fnames{iF});
            this_dataR.(fnames{iF}) = temp(keep{1});
        end
        
        % get the data of interest
        nL = nansum(this_dataR.shuf_perc_odds >= 1-cfg.cutoff); % odds are log(pL/pR) so left means big numbers vs. shuffle
        nR = nansum(this_dataR.shuf_perc_odds <= cfg.cutoff); % odds are log(pL/pR) so right means small numbers vs. shuffle
        nE = length(this_dataR.tvec);
              
        fracL_evt = nL/(nL+nR); fracR_evt = nR/(nL+nR);  % fraction of events (signif events only)
        
        median_z = nanmedian(this_dataR.shuf_z_odds); % median z-score
        median_perc = nanmedian(this_dataR.shuf_perc_odds);  % median percentile vs. shuffle
        
        % append to big data structure
        for iV = 1:length(vars)
            fprintf('..%s\n',vars{iV});
            
            append_this = eval(vars{iV});
            data.(what{iW}).all.(vars{iV}) = cat(1,data.(what{iW}).all.(vars{iV}),append_this); % all rats first
            data.(what{iW}).(this_rat).(vars{iV}) = cat(1,data.(what{iW}).(this_rat).(vars{iV}),append_this);
                     
        end % of variables loop
              
    end % of what loop (all, pre, task, post, etc)
    
end

%% compute means for plotting
nShuf = 1000;
for iR = 1:length(rats)
    for iW = 1:length(what)
        for iV = 1:length(vars)-1 % requires that type is last
        
            this_data = data.(what{iW}).(rats{iR}).(vars{iV});
            this_type = data.(what{iW}).(rats{iR}).this_type;
            
            this_food = this_data(this_type == 0);
            this_water = this_data(this_type == 1);
            
            food_varname = cat(2,'food_',vars{iV});
            water_varname = cat(2,'water_',vars{iV});
            
            data.(what{iW}).(rats{iR}).(food_varname) = this_food;
            data.(what{iW}).(rats{iR}).(water_varname) = this_water;
        
            % shuffles -- randomly reassign food/water for all sessions
            all_food = randi(2, nShuf, 1) - 1;
            shuf_food = []; shuf_water = []; shuf_diff = [];
            for iS = nShuf:-1:1
                
                this_type = this_type(randperm(length(this_type)));
                
                this_food = this_data(this_type == 0); 
                this_water = this_data(this_type == 1); 
                
                shuf_food(iS) = nanmean(this_food);
                shuf_water(iS) = nanmean(this_water);
                shuf_diff(iS) = nanmean(this_food) - nanmean(this_water);
                
            end
            
            foodS_varname = cat(2,'foodShuf_',vars{iV}); foodSsd_varname = cat(2,'foodShufsd_',vars{iV});
            waterS_varname = cat(2,'waterShuf_',vars{iV}); waterSsd_varname = cat(2,'waterShufsd_',vars{iV});
            diffS_varname = cat(2,'diffShuf_',vars{iV}); diffSsd_varname = cat(2,'diffShufsd_',vars{iV}); % food minus water
            
            data.(what{iW}).(rats{iR}).(foodS_varname) = nanmean(shuf_food);
            data.(what{iW}).(rats{iR}).(foodSsd_varname) = nanstd(shuf_food);
            
            data.(what{iW}).(rats{iR}).(waterS_varname) = nanmean(shuf_water);
            data.(what{iW}).(rats{iR}).(waterSsd_varname) = nanstd(shuf_water);
            
            data.(what{iW}).(rats{iR}).(diffS_varname) = nanmean(shuf_diff);
            data.(what{iW}).(rats{iR}).(diffSsd_varname) = nanstd(shuf_diff);
            
        end % of vars
    end
end

%% For splitting vs. replay bias
cfg = []; cfg.sess = [2:6 10:11 13:24];
for i = 1:length(cfg.sess)
    sess_i = cfg.sess(i);
    data.all.all.n_neurons(i) = out{sess_i}.n_neurons;
    data.all.all.FR_diff(i) = mean(out{sess_i}.FR_diff, 'omitnan');
%     data.all.all.FR_diff_zscore(i) = median(abs(out{sess_i}.FR_diff_zscore), 'omitnan');
end

% %% perform PCA and compute representational distance
% num_pcs = 10;
% for i = 1:length(cfg.sess)
%     sess_i = cfg.sess(i);
%     pca_input = out{sess_i}.tc.tc;
%     pca_mean = nanmean(pca_input, 2);
%     pca_input_centered = pca_input - pca_mean;
%     pca_input_centered(isnan(pca_input_centered)) = 0;
% 
%     % PCA
%     [U, S, V] = svd(pca_input_centered);
%     [eigvecs] = U(:, 1:num_pcs);
%     pca_projected = eigvecs' * pca_input_centered;
% 
%     % Compute L vs. R distance in embedding space
%     cp_bin = out{sess_i}.cp_bin;
%     div = ceil(size(pca_projected, 2)/2);
%     pca_L = pca_projected(:, 1:cp_bin);
%     pca_R = pca_projected(:, div:div+cp_bin-1);
%     data.all.all.rep_dist(i) = sum(sqrt((pca_L - pca_R).^2), 'all') / numel(pca_L);
% end

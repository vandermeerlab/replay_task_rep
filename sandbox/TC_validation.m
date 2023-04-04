%% Plot tuning curves
figure;
imagesc(TC.combined.tc.tc); colorbar;
xlabel('bin');
ylabel('neurons');

%% Plot some example trajectories
pos = LoadPos([]);
pos_X = getd(pos,'x');
pos_Y = getd(pos,'y');

figure;
plot(pos_X, pos_Y,'.','Color',[0.7 0.7 0.7],'MarkerSize',4); xlabel('x data'); ylabel('y data');
view(90,90);
hold on;

%% Plot animal's positions where a given neurons fires each spike.
cfg_spikes = {};
cfg_spikes.load_questionable_cells = 1;
S = LoadSpikes(cfg_spikes);
S = restrict(S, metadata.taskvars.trial_iv);
S = restrict(S, run_iv);

iC = 37;
spk_x = interp1(pos.tvec,getd(pos,'x'),S.t{iC},'linear');
spk_y = interp1(pos.tvec,getd(pos,'y'),S.t{iC},'linear');

h = plot(spk_x,spk_y,'.r');

%% Plot position occupancy
figure;
plot(TC.left.tc.occ_hist / nansum(TC.left.tc.occ_hist) * 100);
hold on;
plot(TC.right.tc.occ_hist / nansum(TC.right.tc.occ_hist) * 100);
hold on;
legend('left','right');
xlabel('bin');
ylabel('normalized occupancy');
set(gca,'FontSize',12)
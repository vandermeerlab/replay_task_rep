%% Plot tuning curves
figure;
imagesc(TC.combined.tc.tc); colorbar;
caxis([0, 50]);
xlabel('bin');
ylabel('neurons');

%% Plot some example trajectories
pos = LoadPos([]);
pos_X = getd(pos,'x');
pos_Y = getd(pos,'y');
rot_angle = 90;

figure;
plot(pos_X, pos_Y,'.','Color',[0.7 0.7 0.7],'MarkerSize',4); xlabel('x data'); ylabel('y data');
view(rot_angle,rot_angle);
hold on;

%% Plot animal's positions where a given neurons fires each spike.
cfg_spikes = {};
cfg_spikes.load_questionable_cells = 1;
S = LoadSpikes(cfg_spikes);
S = restrict(S, metadata.taskvars.trial_iv);
S = restrict(S, run_iv);

iC = 62;
spk_x = interp1(pos.tvec,getd(pos,'x'),S.t{iC},'linear');
spk_y = interp1(pos.tvec,getd(pos,'y'),S.t{iC},'linear');

h = plot(spk_x,spk_y,'.r');

%% Plot position occupancy
n_bins = length(unique(TC.left.tc.usr.pos_idx));
left_occupancy = zeros(1, n_bins);
right_occupancy = zeros(1, n_bins);
for i = 1:n_bins
    left_occupancy(i) = sum(TC.left.tc.usr.pos_idx == i);
    right_occupancy(i) = sum(TC.right.tc.usr.pos_idx == i);
end

figure;
plot(left_occupancy / sum(left_occupancy));
hold on;
plot(right_occupancy / sum(right_occupancy));
hold on;
legend('left','right');
xlabel('bin');
ylabel('normalized occupancy');
set(gca,'FontSize',12)
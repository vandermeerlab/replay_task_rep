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

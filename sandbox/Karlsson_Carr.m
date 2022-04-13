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
% center: 1, left: 2, right: 4
pos_trial = zeros(1, length(X));

for pos_i = 1:length(pos_trial)
    % Check if position is in the center arm reward area
    if (X(pos_i) > 60 && X(pos_i) < 90) && (Y(pos_i) > 123 && Y(pos_i) < 169)
        pos_trial(pos_i) = 1;
    % Check if position is in the left arm reward area
    elseif (X(pos_i) > 28 && X(pos_i) < 54) && (Y(pos_i) > 123 && Y(pos_i) < 169)
        pos_trial(pos_i) = 2;
    % Check if position is in the right arm reward area
    elseif (X(pos_i) > 93 && X(pos_i) < 130) && (Y(pos_i) > 123 && Y(pos_i) < 169)
        pos_trial(pos_i) = 4;
    else
        pos_trial(pos_i) = 0;
    end
end

diff_pos_trial = diff(pos_trial);
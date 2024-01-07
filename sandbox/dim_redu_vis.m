%%
food_idx = ~logical(data.all.all.this_type);
water_idx = logical(data.all.all.this_type);

replay_bias_scores = zeros(length(data.all.all.median_z), 1);
replay_bias_scores(food_idx) = -data.all.all.median_z(food_idx);
replay_bias_scores(water_idx) = data.all.all.median_z(water_idx);

[sorted_replay_bias_scores, sorted_idx] = sort(replay_bias_scores);

%%
data = TC;

% Project [L, R] to PCA space.
NumComponents = 10;
for p_i = 1:length(data)
    [proj_data{p_i}, eigvecs{p_i}, pca_mean{p_i}] = perform_pca(data{p_i}, NumComponents);
end

%%
for c_i = 1:ceil(length(sorted_idx) / 6)
    figure;
    for i = 1:6
        s_i = (c_i-1)*6 + i;
        if s_i <= length(sorted_idx)
            subplot(2, 3, i);
            sess_i = sorted_idx(s_i);
            plot_L = plot_3d_trajectory(proj_data{sess_i}.left);
            plot_L.Color = 'r';
            hold on;
            plot_R = plot_3d_trajectory(proj_data{sess_i}.right);
            plot_R.Color = 'b';
            title(sprintf('replay bias: %.2f', sorted_replay_bias_scores(s_i)));
        end
    end
end
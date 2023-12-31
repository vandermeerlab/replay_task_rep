%%
food_idx = ~logical(data.all.all.this_type);
water_idx = logical(data.all.all.this_type);

replay_bias_scores = zeros(length(data.all.all.median_z), 1);
replay_bias_scores(food_idx) = -data.all.all.median_z(food_idx);
replay_bias_scores(water_idx) = data.all.all.median_z(water_idx);

[sorted_replay_bias_scores, sorted_idx] = sort(replay_bias_scores);

%%
TC_matrix = [TC{1}.left, TC{1}.right]';
% labels = [zeros(size(TC{1}.left, 2), 1); ones(size(TC{1}.right, 2), 1)];
TC_umap = run_umap(TC_matrix,'min_dist',0.6,'n_neighbors',50,'metric','cosine','n_components',3);
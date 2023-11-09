%%
n_neurons_per_sess = cellfun(@(x) size(x.left, 1), Q);
min_n_neurons = min(n_neurons_per_sess);

%%
num_subsamples = 100;
num_pcs = 20;
cum_vars = zeros(length(Q), num_pcs, num_subsamples);
min_comps_to_explain = zeros(length(Q), num_subsamples);
for q_i = 1:length(Q)
    for s_i = 1:num_subsamples
        shuffle_indices = randperm(size(Q{q_i}.left, 1));
        % Concatenate Q matrix across left and right trials and perform PCA on it.
        subsample_indices = shuffle_indices(1:min_n_neurons);
        pca_input = [Q{q_i}.left(subsample_indices, :), Q{q_i}.right(subsample_indices, :)];
        pca_mean = mean(pca_input, 2);
        pca_input = pca_input - pca_mean;
        [~, ~, ~, ~, explained, ~] = pca(pca_input);
        cum_vars(q_i, :, s_i) = cumsum(explained(1:num_pcs));
        min_comps_to_explain(q_i, s_i) = find(cum_vars(q_i, :, s_i) > 95, 1);
    end
end

m_comps_to_explain = median(min_comps_to_explain, 2);

%%
figure;
scatter(m_comps_to_explain, abs(data.all.all.median_z), 'filled');

xlabel('Components needed to explain 95% of variance')
ylabel('z-score bias in replay')
title('All data')
set(gca, 'FontSize', 12)

%%
figure;
scatter(m_comps_to_explain(food_idx), data.all.all.median_z(food_idx), 'filled', 'red');
hold on;
scatter(m_comps_to_explain(water_idx), data.all.all.median_z(water_idx), 'filled', 'blue');
xlabel('Components needed to explain 95% of variance')
ylabel('z-score bias in replay')
legend({'food-restricted', 'water-restricted'})
title('All data')
set(gca, 'FontSize', 12)

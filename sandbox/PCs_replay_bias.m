%%
PC_data = TC;
n_neurons_per_sess = cellfun(@(x) size(x.left, 1), PC_data);
min_n_neurons = min(n_neurons_per_sess);

%%
num_pcs = 30;
cum_vars = zeros(length(PC_data), num_pcs);
min_comps_to_explain = zeros(length(PC_data), 1);
non_zero_singulars = zeros(length(PC_data), 1);
part_ratios = zeros(length(PC_data), 1);

for q_i = 1:length(PC_data)
    pca_input = [PC_data{q_i}.left, PC_data{q_i}.right];
    pca_mean = mean(pca_input, 2);
    pca_input = pca_input - pca_mean;
    [coeff,score,latent,tsquared,explained,mu] = pca(pca_input);
%     [~, ~, ~, ~, explained, ~] = pca(pca_input);
    cum_vars(q_i, :) = cumsum(explained(1:num_pcs));
    min_comps_to_explain(q_i) = find(cum_vars(q_i, :) > 95, 1);
    non_zero_singulars(q_i) = sum(latent > 0.01);
    part_ratios(q_i) = sum(latent)^2 / sum(latent.^2);
end

%%
num_subsamples = 100;
num_pcs = 30;
cum_vars = zeros(length(PC_data), num_pcs, num_subsamples);
min_comps_to_explain = zeros(length(PC_data), num_subsamples);
non_zero_singulars = zeros(length(PC_data), num_subsamples);
part_ratios = zeros(length(PC_data), num_subsamples);

for q_i = 1:length(PC_data)
    for s_i = 1:num_subsamples
        shuffle_indices = randperm(size(PC_data{q_i}.left, 1));
        % Concatenate Q matrix across left and right trials and perform PCA on it.
        subsample_indices = shuffle_indices(1:min_n_neurons);
        pca_input = [PC_data{q_i}.left(subsample_indices, :), PC_data{q_i}.right(subsample_indices, :)];
        pca_mean = mean(pca_input, 2);
        pca_input = pca_input - pca_mean;
        [coeff,score,latent,tsquared,explained,mu] = pca(pca_input);
%         [~, ~, ~, ~, explained, ~] = pca(pca_input);
        cum_vars(q_i, :, s_i) = cumsum(explained(1:num_pcs));
        min_comps_to_explain(q_i, s_i) = find(cum_vars(q_i, :, s_i) > 95, 1);
        non_zero_singulars(q_i, s_i) = sum(latent > 0.01);
        part_ratios(q_i, s_i) = sum(latent)^2 / sum(latent.^2);
    end
end

m_comps_to_explain = median(min_comps_to_explain, 2);
m_non_zero_singulars = median(non_zero_singulars, 2);
m_part_ratios = median(part_ratios, 2);

%%
figure;
scatter(m_part_ratios, n_neurons_per_sess, 'filled');

xlabel('Components needed to explain 95% of variance')
ylabel('z-score bias in replay')
title('All data')
set(gca, 'FontSize', 12)

%%
dim_estimates = {m_comps_to_explain, m_non_zero_singulars, m_part_ratios};
method_labels = {'explain 95%', 'non-zero singular values', 'participant ratio'};

for a_i = 1:length(dim_estimates)
    figure;
    dim_est = dim_estimates{a_i};
    scatter(dim_est(food_idx), -data.all.all.median_z(food_idx), 'filled', 'red');
    hold on;
    scatter(dim_est(water_idx), data.all.all.median_z(water_idx), 'filled', 'blue');
    xlabel(method_labels{a_i})
%     xlabel('Components needed to explain 95% of variance')
    ylabel('z-score bias in replay')
%     ylabel('Number of recorded neurons per session')
    legend({'food-restricted', 'water-restricted'})
    set(gca, 'FontSize', 12)
end

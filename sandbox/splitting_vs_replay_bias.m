food_idx = ~logical(data.all.all.this_type);
water_idx = logical(data.all.all.this_type);

%%
scatter(avg_preCP_correct(food_idx), data.all.all.median_z(food_idx), 'filled', 'red');
hold on;
scatter(avg_preCP_correct(water_idx), data.all.all.median_z(water_idx), 'filled', 'blue');
xlabel('Decoding accuracy before choice points (L and R) averaged')
ylabel('z-score bias in replay')
legend({'food-restricted', 'water-restricted'})
title('All data')
set(gca, 'FontSize', 12)

%%
% main plots
epochs = {'pre', 'task', 'post'};
for e_i = 1:length(epochs) % loop over epochs
    figure;
    scatter(avg_preCP_correct(food_idx), data.(epochs{e_i}).all.median_z(food_idx), 'filled', 'red');
    hold on;
    scatter(avg_preCP_correct(water_idx), data.(epochs{e_i}).all.median_z(water_idx), 'filled', 'blue');
    xlabel('Decoding accuracy before choice points (L and R) averaged')
    ylabel('z-score bias in replay')
    legend({'food-restricted', 'water-restricted'})
    title(epochs{e_i})
    set(gca, 'FontSize', 12)
end
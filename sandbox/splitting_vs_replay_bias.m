food_idx = ~logical(data.all.all.this_type);
water_idx = logical(data.all.all.this_type);

% avg_preCP_correct = (nanmean(left_preCP_correct, 2) + nanmean(right_preCP_correct, 2)) / 2;
% avg_postCP_correct = (nanmean(left_postCP_correct, 2) + nanmean(right_postCP_correct, 2)) / 2;

avg_preCP_correct = (left_preCP_correct + right_preCP_correct) / 2;
avg_postCP_correct = (left_postCP_correct + right_postCP_correct) / 2;

avg_entire_correct = (avg_preCP_correct + avg_postCP_correct) / 2;

%%
scatter(avg_correct(food_idx), abs(data.all.all.median_z(food_idx)), 'filled', 'red');
hold on;
scatter(avg_correct(water_idx), abs(data.all.all.median_z(water_idx)), 'filled', 'blue');
xlabel('Decoding accuracy before choice points (L and R) averaged')
% xlabel('Number of recorded neurons per session')
ylabel('z-score bias in replay')
legend({'food-restricted', 'water-restricted'})
title('All data')
set(gca, 'FontSize', 12)

%%
avg_correct_alls = {avg_preCP_correct, avg_postCP_correct, avg_entire_correct};
accuracy_titles = {'pre-CP', 'post-CP', 'entire'};
for a_i = 1:length(avg_correct_alls)
    figure;
    avg_correct = avg_correct_alls{a_i};
    scatter(avg_correct(food_idx), -data.all.all.median_z(food_idx), 'filled', 'red');
    hold on;
    scatter(avg_correct(water_idx), data.all.all.median_z(water_idx), 'filled', 'blue');
    xlabel('Decoding accuracy')
%     ylabel('Number of recorded neurons per session')
    ylabel('opposite-to-behavior bias in replay')
    legend({'food-restricted', 'water-restricted'})
%     xlim([0.4, 1])
    title(accuracy_titles{a_i})
    set(gca, 'FontSize', 14)
    save_as_eps('~/Desktop', [accuracy_titles{a_i}, ''])
end
%%
avg_correct_alls = {avg_preCP_correct, avg_postCP_correct, avg_entire_correct};
epochs = {'pre', 'task', 'post'};
figure;
for a_i = 1:length(avg_correct_alls)
    avg_correct = avg_correct_alls{a_i};
    for e_i = 1:length(epochs)
        subplot(3, 3, (a_i-1)*3+e_i)
        scatter(avg_correct(food_idx), -data.(epochs{e_i}).all.median_z(food_idx), 'filled', 'red');
        hold on;
        scatter(avg_correct(water_idx), data.(epochs{e_i}).all.median_z(water_idx), 'filled', 'blue');
        xlabel('Decoding accuracy')
        ylabel('z-score bias in replay')
        legend({'food-restricted', 'water-restricted'})
        title(epochs{e_i})
        set(gca, 'FontSize', 12)
    end
end
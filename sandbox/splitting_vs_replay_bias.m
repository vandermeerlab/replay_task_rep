food_idx = ~logical(data.all.all.this_type);
water_idx = logical(data.all.all.this_type);

replay_bias_scores = zeros(length(data.all.all.median_z), 1);
replay_bias_scores(food_idx) = -data.all.all.median_z(food_idx);
replay_bias_scores(water_idx) = data.all.all.median_z(water_idx);

% avg_preCP_correct = (nanmean(left_preCP_correct, 2) + nanmean(right_preCP_correct, 2)) / 2;
% avg_postCP_correct = (nanmean(left_postCP_correct, 2) + nanmean(right_postCP_correct, 2)) / 2;

avg_preCP_correct = (left_preCP_correct + right_preCP_correct) / 2;
avg_postCP_correct = (left_postCP_correct + right_postCP_correct) / 2;

avg_entire_correct = (avg_preCP_correct + avg_postCP_correct) / 2;

%%
figure;
p = polyfit(avg_preCP_correct, replay_bias_scores, 1);
predicted_bias = polyval(p, avg_preCP_correct);

scatter(avg_preCP_correct, replay_bias_scores, 'filled');
hold on;
plot(avg_preCP_correct, predicted_bias, '--', 'Color', [.8 .8 .8], 'LineWidth', 2); 

xlabel('Decoding accuracy (L vs. R)')
ylabel('opposite-to-behavior bias in replay')
title(['r = ', num2str(R(2, 1))])
set(gca, 'FontSize', 14)

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
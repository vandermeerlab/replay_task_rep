%% Ensemble splitting
left_preCP_correct = cellfun(@(x) x.left_preCP_correct, out);
left_postCP_correct = cellfun(@(x) x.left_postCP_correct, out);
right_preCP_correct = cellfun(@(x) x.right_preCP_correct, out);
right_postCP_correct = cellfun(@(x) x.right_postCP_correct, out);

left_chance_pre = cellfun(@(x) x.left_chance_pre, out);
left_chance_post = cellfun(@(x) x.left_chance_post, out);
right_chance_pre = cellfun(@(x) x.right_chance_pre, out);
right_chance_post = cellfun(@(x) x.right_chance_post, out);

ywhats = {left_preCP_correct, left_postCP_correct, right_preCP_correct, right_postCP_correct};
ychances = {left_chance_pre, left_chance_post, right_chance_pre, right_chance_post};
cols = 'rrbb';
xd = [1 2 4 5];

figure;

for iW = 1:length(ywhats)
    
   ydata = ywhats{iW}; ychance = ychances{iW};
   h = bar(xd(iW),nanmean(ydata)); set(h,'FaceColor',cols(iW),'EdgeColor','none','BarWidth',0.6);
   hold on;
   h2 = plot(repmat(xd(iW),[1 length(ydata)]),ydata,'.k');
   h3 = plot([xd(iW)-0.4 xd(iW)+0.4],[nanmean(ychance) nanmean(ychance)],':k','LineWidth',2);
    
   % report means and SEMs; statistics vs. chance
   m = nanmean(ydata); s = nanstd(ydata)./sqrt(4); p = ranksum(ydata,ychance);
   fprintf('%s: %.2f +/- %.2f (SEM), p = %2.2e vs. chance\n',ywhats{iW},m,s,p);
   
end
set(gca,'LineWidth',1,'TickDir','out','FontSize',18,'XLim',[0 6],'XTick',xd,'XTickLabel',{'preL','postL','preR','postR'},'YLim',[0 1],'YTick',0:0.25:1);
box off;

%% Single-cell splitting
FR_diff = cellfun(@(x) mean(x.left_preCP_correct, 'omitnan'), out);

%% Splitting vs. replay bias
avg_preCP_correct = (left_preCP_correct + right_preCP_correct) / 2;

x = avg_preCP_correct(~isnan(oppo_replay_bias));
y = oppo_replay_bias(~isnan(oppo_replay_bias));

figure;
p = polyfit(x, y, 1);
predicted_y = polyval(p, x);

[corr, p_val] = corrcoef(x, y);
scatter(x, y, 'filled', 'MarkerFaceColor', [105/255 105/255 105/255]);
hold on;
plot(x, predicted_y, '--', 'Color', [198/255 113/255 113/255], 'LineWidth', 2); 

xlabel('ensemble splitter strength')
ylabel('$\left|FR_{left} - FR_{right}\right|$ ', 'Interpreter','latex')
title(['r = ', num2str(corr(2, 1)), '; p-value = ', num2str(p_val(2, 1))])
% title(['r = ', num2str(bias_corr(2, 1))])

set(gca, 'FontSize', 24)
data1 = pMetAndAllFDR(1:length(pVals));
data2 = pMetAndAllFDR(length(pVals):end);

colors = brewermap(5,'PRGn');
alpha(.95);

f = figure;
c = cdfplot(log10(data1));
set(c,'Color',colors(end,:),'LineWidth',2.5);
hold on;
c = cdfplot(log10(data2));
set(c,'Color',colors(1,:),'LineWidth',2.5);

ranksum(data1,data2)

lgd = legend('Metabolic genes','Non-metabolic genes');
lgd.FontSize = 18;
legend('Location','west');

xlab = xlabel('log10 p-values DMEM vs RPMI');% x-axis label;
ylab = ylabel('Cumulative frequency'); % y- axis label

%x axis limits
xlim([-6 -1.301]);
ylim([0 0.4]);
%axes font size
set(gca,'FontSize',18);
% title and font size
hTitle = title('');
%set(hTitle,'FontSize',18);
set(gcf,'Position',[100 100 800 600])
set(gca,'FontSize',16);

global Config;
folder = Config.FIGURES_FOLDER;

filePath = fullfile(folder,sprintf('%s.png','1c'));
print(gcf,filePath,'-dpng','-r300');
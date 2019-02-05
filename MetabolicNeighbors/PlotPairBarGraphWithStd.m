function PlotPairBarGraphWithStd(countsFound,meansRand,stdRand,labels,figName)
figure;

b = bar([countsFound,meansRand])

set(gca,'xticklabel',labels)% Sets tick labels back on the Axis
set(gca,'FontSize',16);
xtickangle(30);
ylabel('Fraction of genes correlated to expression');
hold on;
for i=1:length(countsFound)
    e = errorbar(i+0.15,meansRand(i),stdRand(i));
    e.Color = 'black';
end
lgd = legend('Metabolic genes','Random genes');
lgd.FontSize = 16;
lgd.Location = 'northwest';

colors1 = brewermap(5,'PRGn');
colors2 = brewermap(5,'PiYG');

set(b(1),'FaceColor',colors1(end,:),'facealpha',.95,'edgecolor','none');
set(b(2),'FaceColor',colors2(1,:),'facealpha',.95,'edgecolor','none');
ylim([0, 0.75]);
grid on;
set(gcf,'Position',[100 100 800 600])

global Config;
folder = Config.FIGURES_FOLDER;

filePath = fullfile(folder,sprintf('%s.png',figName));
print(gcf,filePath,'-dpng','-r300');
end


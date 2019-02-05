function [ ] = CreateBarGraphMostSensitiveGenes(pVals,pMetAndAllFDR,genesRes,r,fileName)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

pVals = pMetAndAllFDR(1:length(pVals));

ids = find(pVals<0.05);
pValsSelected = pVals(ids);
genesSelected = genesRes(ids);
r = r(ids);

[vals,I] = sort(pValsSelected,'desc');
vals = -1*log10(vals);
genesSelected = genesSelected(I);
r = r(I);

%figure
%b =bar(vals,'FaceColor','flat');

colors = brewermap(5,'RdBu');

mydata=vals;
figure(1)
hold on
for i = 1:length(mydata)
    h=barh(i,mydata(i));
    if r(i) < 0
        set(h,'FaceColor',colors(1,:),'facealpha',.8,'edgecolor','none');
    else
        set(h,'FaceColor',colors(end,:),'facealpha',.8,'edgecolor','none');
    end
    barsArr(i)= h(1);
end
%hold off

set(gca, 'YTick', []);
a = 1:length(genesSelected);
set(gca, 'YTickLabel', genesSelected);
set(gca, 'YTick', a)
%xtickangle(90);
%lgd = legend(titles);
%lgd.FontSize = 30;
xlabel('-log10(p-value)');

set(gcf,'Position',[0 0 500 40*length(genesSelected)+60])
xlim([0,20])
grid on;
set(gca,'fontsize',18);

lgd= legend(barsArr([end-1,end]),'RPMI','DMEM');
lgd.Location = 'southeast';
lgd.FontSize = 18;
global Config;
filePath = fullfile(Config.FIGURES_FOLDER,sprintf('%s.png',fileName));
print(gcf,filePath,'-dpng','-r300');
end


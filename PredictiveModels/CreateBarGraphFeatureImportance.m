function [  ] = CreateBarGraphFeatureImportance(res,modelMet,fileName)
%CreatePieChartFeatureImportance
classes = {'ExpCurrent','ExpRelated','ExpRelated','CnvCurrent','CnvRelated','CnvRelated','Media','Mut','CancerType'};
titles = {'Exp-Self','Exp-Iso','Exp-Related','Cnv-Self','Cnv-Iso','Cnv-Related','Media','Mutation','Lineage'};
sortFeats = 1;
pFDR= res.pValuesFDR;
bestImpFeatures = res.bestImpFeatures;
ids = find(pFDR<0.05);
numTopArr = [1,5];
counters = zeros(length(classes),length(numTopArr));

for k=1:length(numTopArr)
    for i=1:length(ids)
        curId = ids(i);
        curFeatures = bestImpFeatures{curId};
        feat= min(length(curFeatures),numTopArr(k));
        curFeatures = curFeatures(1:feat);
        for j=1:length(classes)
            %Check for Iso
            boolFound = 0;

            if(j==2 || j==3 || j==5 ||j==6)
                targetGene = res.genes{curId};
                for t=1:length(curFeatures)
                    splited = strsplit(curFeatures{t},'_');
                    if(length(splited)>1)
                        featureGene = splited{2};
                        [sameReact,iso,complex] = CheckIfIsoOrComplex(targetGene,featureGene,modelMet);
                        if(sameReact && (j==2 ||j==5))
                            if(contains(curFeatures{t},classes{j}))
                                boolFound = 1;
                            end
                        elseif(~sameReact &&(j==3 ||j==6))
                            if(contains(curFeatures{t},classes{j}))
                                boolFound = 1;
                            end
                        end
                    end
                end
                
            else                
                if(sum(contains(curFeatures,classes{j})))
                    boolFound = 1;
                end
            end
            counters(j,k) = counters(j,k)+boolFound;
        end
    end
end
freqs = counters/length(ids);
figure
if(sortFeats)
    [sorted,I] = sort(freqs(:,1),'descend');
    freqs = freqs(I,:);
    titles = titles(I);
end

b = bar(freqs*100);
set(gca, 'XTickLabel', titles);
set(gca,'fontsize',12);
title('');

grid on;

lgd = legend('Top predictive feature','Top 5 predictive features');
lgd.FontSize = 16;
lgd.Location = 'northeast';
colors = brewermap(5,'GnBu');
set(b(1),'FaceColor',colors(2,:),'facealpha',.99,'edgecolor','none');
set(b(2),'FaceColor',colors(4,:),'facealpha',.99,'edgecolor','none');
xtickangle(30);

yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis

set(gca,'FontSize',11);
ylabel('Percent of models');
set(gcf,'Position',[100 100 800 600])
set(gca,'fontsize',16);

global Config;
folder = Config.FIGURES_FOLDER;
filePath = fullfile(folder,sprintf('%s.png',fileName));
print(gcf,filePath,'-dpng','-r300');

end


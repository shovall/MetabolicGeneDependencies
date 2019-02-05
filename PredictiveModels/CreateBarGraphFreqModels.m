function CreateBarGraphFreqModels(resRelated,resRelatedCancer,depMapMetadata,genes)

if(nargin<4)
    genes = resRelated.genes;
end
[exists,locs] = ismember(genes,resRelated.genes);
locsRelated = locs(exists);
num1FDRAlone = length(find(resRelated.pValuesFDR(locsRelated)<0.05));

[exists,locs] = ismember(genes,resRelated.genes);
locsCancer = locs(exists);
num2FDRAlone = length(find(resRelatedCancer.pValuesFDR(locsCancer)<0.05));

if(isstruct(depMapMetadata))
    [exists,locs] = ismember(genes,depMapMetadata.genes);
    numUnbiased = length(find(depMapMetadata.unbiasedModelExists(locs(exists))));
	numRelated = length(find(depMapMetadata.relatedModelExists(locs(exists))));
    vals = [numUnbiased;numRelated;num1FDRAlone;num2FDRAlone];
    numAnalyzedGenesDepMap = length(find(depMapMetadata.isAnalyzed(locs(exists))));
    vals(1:2)= vals(1:2)/numAnalyzedGenesDepMap*100;
    vals(3:end)= vals(3:end)./length(locsRelated)*100;
else
    vals = [num1FDRAlone;num2FDRAlone];
    vals=vals./length(locsRelated)*100;
end

mydata=vals(1:end);

figure(1)
grid on;
box on;
hold on
for i = 1:length(mydata)
    h=bar(i,mydata(i));
    if (i==1 || i==2) && isstruct(depMapMetadata)
		colors = brewermap(3,'Set2');
        set(h,'FaceColor',colors(2,:),'facealpha',.9,'edgecolor','none');
    else
		colors = brewermap(5,'BrBG');
		if(isstruct(depMapMetadata))
			color = colors(end,:);
		else
			color = colors(1,:);
		end
        set(h,'FaceColor',color,'facealpha',.9,'edgecolor','none');
    end
end

if(isstruct(depMapMetadata))
    a = 1:4;
    titles = {'DepMap Unbiased','DepMap related','Related & Media','Related, Media & Lineage'};
else
    a = 1:2;
    titles = {'Related & Media','Related, Media & Lineage'};
    
end
%ylim([0 33]);

set(gca, 'XTick', []);
set(gca, 'XTickLabel', titles);
set(gca, 'XTick', a)

set(gca,'fontsize',16);
title('');%Fraction of genes with predictive models');
ylabel('Percent of models');
if(isstruct(depMapMetadata))
    set(gcf,'Position',[100 100 800 600])
else
    set(gcf,'Position',[100 100 800*1 600])
end

yticks = [get(gca,'ytick')]'; % There is a transpose operation here.
percentsy = repmat('%', length(yticks),1);  %  equal to the size
yticklabel = [num2str(yticks) percentsy]; % concatenates the tick labels
set(gca,'yticklabel',yticklabel)% Sets tick labels back on the Axis
xtickangle(15);

if(isstruct(depMapMetadata))
    if(nargin<4)
        fileName = '3a';
    else
        fileName = 'SupModels6SigRNAi';
    end
else
    if(nargin<4)
        fileName = '3b';
    else
        fileName = 'SupModels6SigCRISPR';
    end
end

global Config;
folder = Config.FIGURES_FOLDER;
filePath = fullfile(folder,sprintf('%s.png',fileName));
print(gcf,filePath,'-dpng','-r300');

end

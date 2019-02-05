function [genes,pVals,pValsFDR,pValsSpearmanPartialCond,rVals,T] = GetTableMediaScreening(plotHist,analysisTitle,genesToCheck,depMapData,depMapMetadata,mediaTypesCells,mediaTypesSelected,thresh)
titlesForTable = {};
mediaInds = cell(1,length(mediaTypesSelected));
mediaIndsDict = zeros(length(mediaTypesSelected));
conditions = cell(1,length(mediaTypesSelected));
cancerTypes = cell(1,length(mediaTypesSelected));

k=1;
for i=1:length(mediaTypesSelected)
    [exists,inds] = ismember(upper(mediaTypesCells{i}), upper(depMapData.celllines));
    inds = inds(exists);
    mediaInds{i} = inds;
    [exists,inds] = ismember(upper(mediaTypesCells{i}), upper(depMapMetadata.celllines));
    inds = inds(exists);
    if(~isempty(depMapMetadata.conditions))
        conditions{i} = depMapMetadata.conditions(inds);
    end
    if(~isempty(depMapMetadata.primaryDisease))
        cancerTypes{i} = depMapMetadata.primaryDisease(inds);
    end
    
    
    for j=i+1:length(mediaTypesSelected)
        titlesForTable{end+1} = strcat(mediaTypesSelected{i},'vs',mediaTypesSelected{j});
        mediaIndsDict(i,j) = k;
        k= k+1;
    end
end

if(~isempty(cancerTypes{1}))
    cancerTypesVec = [cancerTypes{1}',cancerTypes{2}'];
    cancerTypesVec = cellfun(@(x) strip(x),cancerTypesVec,'un',0);
    [cancerTypesVecGrouped,cancersType] = findgroups(cancerTypesVec);
end


if(~isempty(conditions{1}))
    condVec = [conditions{1}',conditions{2}'];
    condVecGrouped = findgroups(condVec);
    %adherent is 1
    [exists,~] = ismember(strtrim(upper(condVec)),'SUSPENSION');
    condVecGrouped(find(exists))=0;
    condVecGrouped(find(condVecGrouped>1))=NaN;
end

[ genesEssen,idsInDepMap,amountsCellLines ] = FindGenesSigmaThresh( depMapData, thresh,0 );
[isExists,locs] = ismember(genesToCheck,genesEssen);
idsInDepMapSelected = idsInDepMap(locs(isExists));

pVals = zeros(length(idsInDepMapSelected),length(titlesForTable))+NaN;
rVals = zeros(length(idsInDepMapSelected),length(titlesForTable));
zVals = zeros(length(idsInDepMapSelected),length(titlesForTable));

pValsSpearman = zeros(length(idsInDepMapSelected),length(titlesForTable))+NaN;
rSpearman= zeros(length(idsInDepMapSelected),length(titlesForTable));

pValsSpearmanPartialCond = zeros(length(idsInDepMapSelected),length(titlesForTable))+NaN;
rSpearmanPartialCond= zeros(length(idsInDepMapSelected),length(titlesForTable));
pValsSpearmanPartialCancer = zeros(length(idsInDepMapSelected),length(titlesForTable))+NaN;
rSpearmanPartialCancer= zeros(length(idsInDepMapSelected),length(titlesForTable));


for k=1:length(idsInDepMapSelected)
    index = idsInDepMapSelected(k);
    col = 1;
    for i=1:length(mediaTypesSelected)
        data1 = depMapData.data(index,mediaInds{i});
        for j=i+1:length(mediaTypesSelected)
            data2 = depMapData.data(index,mediaInds{j});
            if(~(sum(~isnan(data1))==0 || sum(~isnan(data2))==0))
                [p,~,stat] = ranksum(data1, data2);
                pVals(k,col) = p;
                rVals(k,col) = median(data1)-median(data2);
                zVals(k,col) = stat.zval;
                
                concData =[data1,data2]';
                isNotNanData1 = ~isnan(data1);
                isNotNanData2 = ~isnan(data2);
                mediaVec = [zeros(1,sum(isNotNanData1)),zeros(1,sum(isNotNanData2))+1]';
                isNotNan = ~isnan(concData);
                concData = concData(isNotNan);
                
                [rho,p] = corr(concData,mediaVec,'type','Spearman');
                pValsSpearman(k,col) = p;
                rSpearman(k,col) = rho;
                
                if(~isempty(depMapMetadata.conditions))
                condVecAfterNans = condVecGrouped(isNotNan)';               
                isNotNanCond = ~isnan(condVecAfterNans);          
                [rho,p] = partialcorr(concData(isNotNanCond),mediaVec(isNotNanCond),condVecAfterNans(isNotNanCond),'type','Spearman');
                pValsSpearmanPartialCond(k,col) = p;
                rSpearmanPartialCond(k,col) = rho;
                end
                
                if(~isempty(depMapMetadata.primaryDisease))
                cancerVecAfterNans = cancerTypesVecGrouped(isNotNan)';
                [rho,p] = partialcorr(concData,mediaVec,cancerVecAfterNans,'type','Spearman');
                pValsSpearmanPartialCancer(k,col) = p;
                rSpearmanPartialCancer(k,col) = rho;
                end
            end
            col=col+1;
        end
    end
end
pValsFDR = CalcFDR(pVals);

if(plotHist)
    plotAllHistograms(analysisTitle,pValsFDR,idsInDepMapSelected,depMapData,mediaInds,mediaIndsDict,mediaTypesSelected);
end

idsTable = array2table([1:length(idsInDepMapSelected)]','VariableNames',{'id'});
tabeleGenes = cell2table(depMapData.genes(idsInDepMapSelected),'VariableNames',{'Gene'});
[tableLinksToHist] = GetTableLinks(pValsFDR,analysisTitle);
pValTable = array2table(pVals,'VariableNames',titlesForTable);
pValFDRTable = array2table(pValsFDR,'VariableNames',strcat(titlesForTable,'_FDR'));
rValTable = array2table([rVals,zVals],'VariableNames',{'r','z'});
pSpear = array2table(pValsSpearman,'VariableNames',{'pSpearman'});
pSpearFDR = array2table(CalcFDR(pValsSpearman),'VariableNames',{'pSpearmanFDR'});
rSpear = array2table(rSpearman,'VariableNames',{'rSpearman'});
pSpearPartCond = array2table(pValsSpearmanPartialCond,'VariableNames',{'pSpearmanPartialCond'});
pSpearPartCondFDR = array2table(CalcFDR(pValsSpearmanPartialCond),'VariableNames',{'pSpearmanPartialCondFDR'});
rSpearPartcond= array2table(rSpearmanPartialCond,'VariableNames',{'rSpearmanPartialCond'});
pSpearPartCancer = array2table(pValsSpearmanPartialCancer,'VariableNames',{'pSpearmanPartialCancer'});
pSpearPartCancerFDR = array2table(CalcFDR(pValsSpearmanPartialCancer),'VariableNames',{'pSpearmanPartialCancerFDR'});
rSpearPartCancer= array2table(rSpearmanPartialCancer,'VariableNames',{'rSpearmanPartialCancer'});

T = [idsTable,tabeleGenes,tableLinksToHist,pValTable,pValFDRTable,rValTable,pSpear,pSpearFDR,rSpear,pSpearPartCond,pSpearPartCondFDR,rSpearPartcond,pSpearPartCancer,pSpearPartCancerFDR,rSpearPartCancer];
genes = depMapData.genes(idsInDepMapSelected);
end
function [tableLinksToHist] = GetTableLinks(pValsFDR,analysisTitle)
histLinks = cell(length(pValsFDR),1);

for i=1:size(pValsFDR,1)
    ids = find(pValsFDR(i)<0.05);
    for j=ids'
        plotFolderName = 'Histograms';
        filePath = sprintf('%s/%s-%d-%d.png',plotFolderName,analysisTitle,i,j);
        histLinks{i} = sprintf('=HYPERLINK("%s","hist")',filePath);
    end
end
tableLinksToHist = cell2table(histLinks,'VariableNames',{'histogram'});
end
function [] = plotAllHistograms(analysisTitle,pValsFDR,idsInDepMapSelected,depMapData,mediaInds,mediaIndsDict,mediaTypesSelected)

for k=8:size(pValsFDR,1)
    idsSignif = find(pValsFDR(k)<0.05);
    index = idsInDepMapSelected(k);
    for i=idsSignif'
        [t1,t2] = ind2sub(size(mediaIndsDict),find(mediaIndsDict==i));
        
        data1 = depMapData.data(index,mediaInds{t1});
        data2 = depMapData.data(index,mediaInds{t2});
        plotHistogram(data1,data2,analysisTitle,k,i,mediaTypesSelected{t1},mediaTypesSelected{t2});
    end
end

end
function [] = plotHistogram(data1,data2,analysisTitle,fileId,fileId2,type1,type2)
global Config
colors = brewermap(5,'RdBu');
f = figure('Visible','Off');
hold off;
h1 = histogram(data1, 'facecolor',colors(1,:),'facealpha',.7,'edgecolor','none','Normalization','probability');
hold on;
h2 = histogram(data2, 'facecolor',colors(end,:),'facealpha',.7,'edgecolor','none','Normalization','probability');
h1.BinWidth = h2.BinWidth;
xlabel('Dependency score');
ylabel('Fraction of cell lines');
ylim([0 0.7])
lgd = legend([h2(1) h1(1)], type2, type1);
lgd.FontSize = 18;
lgd.Location = 'northwest';
%axes font size
set(gcf,'Position',[100 100 800 600])
set(gca,'FontSize',18);
grid on
path = sprintf('%s\\%s-%d-%d.png',Config.MEDIA_SCREENING_HISTOGRAMS,analysisTitle,fileId,fileId2);
print(gcf,path,'-dpng','-r600');
end

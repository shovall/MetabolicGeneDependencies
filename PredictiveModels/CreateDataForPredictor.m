function [ res ] = CreateDataForPredictor(gene, profileDataCCLE,ccleCNV, driverMutCCLEData,depMapData,depMapMetadata,takeAllGenes,resDist,randomNeighbors)
if(nargin<9)
    randomNeighbors = 0;
end
global Config;

maxCurrencyThresh = Config.CURRENCY_MET_THRESH; 
numGenesRequired = Config.MAX_NUM_OF_NEIGHBORS;
nansFeatureThresh = Config.FRAC_NAN_IGNORE_FEATURE;

[~,locGeneDepMap] = ismember(gene,depMapData.genes);
idsCellsDepMap = ~isnan(depMapData.data(locGeneDepMap,:));
cellNamesSelected = depMapData.celllines(idsCellsDepMap);
idsCellsDepMap = find(idsCellsDepMap);

[exists,locsMutData] = ismember(cellNamesSelected,driverMutCCLEData.caseId);
locsMutData = locsMutData(exists);
cellNamesSelected = cellNamesSelected(exists);
idsCellsDepMap = idsCellsDepMap(exists);

[exists,locsCellsExpData] = ismember(cellNamesSelected,upper(profileDataCCLE.caseId));
locsMutData = locsMutData(exists);
locsCellsExpData = locsCellsExpData(exists);
cellNamesSelected = cellNamesSelected(exists);
idsCellsDepMap = idsCellsDepMap(exists);

[exists,locsCellsCnvData] = ismember(cellNamesSelected,upper(ccleCNV.caseId));
locsMutData = locsMutData(exists);
locsCellsExpData = locsCellsExpData(exists);
cellNamesSelected = cellNamesSelected(exists);
idsCellsDepMap = idsCellsDepMap(exists);
locsCellsCnvData = locsCellsCnvData(exists);

%Mutations Data
driverGenesData = driverMutCCLEData.discreteData(:,locsMutData)';

% Essentiality vector
res.Z = depMapData.data(locGeneDepMap,idsCellsDepMap)';

[~,locsMetadata] = ismember(cellNamesSelected,depMapMetadata.celllines);
% Media and cancer type data
mediaTypeData = depMapMetadata.mediaTypeGrouped(locsMetadata);
[cancerTypeDataGrouped,cancerTypeNames] = findgroups(strtrim(depMapMetadata.primaryDisease(locsMetadata)));

selectedCancers = [];
cancerTypeData = [];
for i=1:max(cancerTypeDataGrouped)
    ids = find(cancerTypeDataGrouped==i);
    if(length(ids)>=Config.CANCER_LINEAGE_MIN_NUMBER) 
        selectedCancers(end+1) = i;
        curCancerData = zeros(length(cancerTypeDataGrouped),1);
        curCancerData(ids) = 1;
        cancerTypeData = [cancerTypeData,curCancerData];
    end
end
cancerTypes = cancerTypeNames(selectedCancers);

if(~takeAllGenes)
    [~,loc] = ismember(gene,resDist.genes);
    [sorted,~] = sort(resDist.data{loc}.dist);
    thresh = max(sorted(1:numGenesRequired));
    thresh  = min(thresh,maxCurrencyThresh);
    locRelatedGenesThresh = find(resDist.data{loc}.dist<=thresh);
    relatedGenes = resDist.data{loc}.genes(locRelatedGenesThresh);    
    if(randomNeighbors)
        [~,locsGenesInExp] = ismember(relatedGenes,profileDataCCLE.common);
        numGenesUsed = length(find(locsGenesInExp));
        randInds = randperm(length(profileDataCCLE.common),numGenesUsed);
        relatedGenes = profileDataCCLE.common(randInds);
    end
else
    relatedGenes = union(profileDataCCLE.common,ccleCNV.common);
    relatedGenes = setdiff(relatedGenes,gene);
end
[~,locsGenesInExp] = ismember(relatedGenes,profileDataCCLE.common);
[~,locCurrentGene] = ismember(gene,profileDataCCLE.common);

locsGenesInExp = [locCurrentGene;locsGenesInExp];
idsUsed = [];
for i=1:length(locsGenesInExp)
    if(locsGenesInExp(i)>0)
        curData = profileDataCCLE.data(locsGenesInExp(i),locsCellsExpData);
        if(~(sum(isnan(curData))>nansFeatureThresh*length(curData)))
            idsUsed = [idsUsed,i];
        end
    end
end
genesIdsUsed = locsGenesInExp(idsUsed);
varNames = [];
if(~isempty(idsUsed) && idsUsed(1) == 1)
    varNames = strcat('ExpCurrent_',profileDataCCLE.common(locsGenesInExp(1)));
    varNames = [varNames,strcat('ExpRelated_',profileDataCCLE.common(genesIdsUsed(2:end)))'];
else
    varNames = [varNames,strcat('ExpRelated_',profileDataCCLE.common(genesIdsUsed))'];
end
expData = profileDataCCLE.data(locsGenesInExp(idsUsed),locsCellsExpData)';

[exists,locsGenesInCnv] = ismember(relatedGenes,ccleCNV.common);
locsGenesInCnv = locsGenesInCnv(exists);
[~,locCurrentGene] = ismember(gene,ccleCNV.common);

locsGenesInCnv = [locCurrentGene;locsGenesInCnv];
idsUsed = [];
for i=1:length(locsGenesInCnv)
    if(locsGenesInCnv(i)>0)
        curData = ccleCNV.data(locsGenesInCnv(i),locsCellsCnvData);
        if(~(sum(isnan(curData))>nansFeatureThresh*length(curData)))
            idsUsed = [idsUsed,i];
        end
    end
end
genesIdsUsed = locsGenesInCnv(idsUsed);
if(~isempty(idsUsed) && idsUsed(1) == 1)
    varNames = [varNames,strcat('CnvCurrent_',ccleCNV.common(locsGenesInCnv(1)))];
    varNames = [varNames,strcat('CnvRelated_',ccleCNV.common(genesIdsUsed(2:end)))'];
else
    varNames = [varNames,strcat('CnvRelated_',ccleCNV.common(genesIdsUsed))'];
end
cnvData = ccleCNV.data(locsGenesInCnv(idsUsed),locsCellsCnvData)';

% Including CNV of oncogenes
[~,locsGenesInCnv] = ismember(driverMutCCLEData.common,ccleCNV.common);
idsUsed = [];
for i=1:length(locsGenesInCnv)
    if(locsGenesInCnv(i)>0)
        curData = ccleCNV.data(locsGenesInCnv(i),locsCellsCnvData);
        if(~(sum(isnan(curData))>nansFeatureThresh*length(curData)))
            idsUsed = [idsUsed,i];
        end
    end
end
oncoCnvGeneLocs = locsGenesInCnv(idsUsed);
oncoCnvData = ccleCNV.data(locsGenesInCnv(idsUsed),locsCellsCnvData)';

res.data = [expData,cnvData,mediaTypeData,cancerTypeData,driverGenesData,oncoCnvData];

varNames = [varNames,'Media'];
varNames = [varNames,strcat('CancerType_',cancerTypes)'];
varNames = [varNames,strcat('Mut_',driverMutCCLEData.common)'];
varNames = [varNames,strcat('OncoCnv_',ccleCNV.common(oncoCnvGeneLocs))'];

res.varNames = varNames;
res.isCategorized = zeros(size(varNames));
res.isCategorized(size(expData,2)+1:end) = 1;
res.idsCellsDepMap = idsCellsDepMap;
res.cellNamesSelected = cellNamesSelected;
end


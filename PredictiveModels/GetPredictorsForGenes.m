function [res] = GetPredictorsForGenes( genesSelected ,profileDataCCLE, ccleCNV,driverMutCCLEData,depMapData,depMapMetadata,varTypes,takeAllGenes,resDist,equalWeights,randomNeighbors)
if(nargin<11)
    randomNeighbors = 0;
end
bestImpFeatures = cell(length(genesSelected),1);
bestImpFeaturesRank = cell(length(genesSelected),1);
pValues = zeros(length(genesSelected),1)+1;
rValues = zeros(length(genesSelected),1);

parfor i=1:length(genesSelected)
    if(mod(i,50)==0)
        disp(i);
    end
    res = CreateDataForPredictor(genesSelected{i}, profileDataCCLE,ccleCNV,...
        driverMutCCLEData,depMapData,depMapMetadata,takeAllGenes,resDist,randomNeighbors);
    if(isempty(res.data))
        continue;
    end

    features = find(contains(res.varNames,varTypes));

    
    if(isempty(features))
        pValues(i) = nan;
        rValues(i) = nan;
        continue;
    end
    
    X = res.data(:,features);
    Y = res.Z;
    varNames = res.varNames(:,features);
    isCategorized = res.isCategorized(features);
        
    [r,p,featuresSelectedIds,imp,EnsembleCV] = GetPredictor(X,Y,isCategorized,~randomNeighbors,equalWeights);    
    if(isnan(imp))
        sortedImp = NaN;
        sortedImpFeat = NaN;
    else 
        [sortedImp,I] = sort(imp,'descend');
        featNamesSorted= varNames(featuresSelectedIds(I));
        numImportantFeat = length(find(sortedImp>0));
        sortedImp = sortedImp(1:numImportantFeat);
        sortedImpFeat = featNamesSorted(1:numImportantFeat);
    end
    bestImpFeatures{i} = sortedImpFeat;
    pValues(i) = p;
    rValues(i) = r;
    bestImpFeaturesRank{i} = sortedImp;
end

res = struct;
res.genes = genesSelected;
res.bestImpFeatures = bestImpFeatures;
res.bestImpFeaturesRank = bestImpFeaturesRank;
res.pValues = pValues;
res.rValues = rValues;
end


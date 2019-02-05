function [res] = GetTableMediaScreeningCancerCancel(genesToCheck,depMapData,depMapMetadata,mediaTypesCells,mediaTypesSelected,thresh)
titlesForTable = {};
mediaInds = cell(1,length(mediaTypesSelected));
mediaIndsDict = zeros(length(mediaTypesSelected));
conditions = cell(1,length(mediaTypesSelected));
cancerTypes = cell(1,length(mediaTypesSelected));

k=1;
for i=1:length(mediaTypesSelected)
    [exists,inds] = ismember(mediaTypesCells{i}, upper(depMapData.celllines));
    inds = inds(exists);
    mediaInds{i} = inds;
    [exists,inds] = ismember(mediaTypesCells{i}, upper(depMapMetadata.celllines));
    inds = inds(exists);
    if(~isempty(depMapMetadata.conditions))
        conditions{i} = depMapMetadata.conditions(inds);
    end
    cancerTypes{i} = depMapMetadata.primaryDisease(inds);
    
    for j=i+1:length(mediaTypesSelected)
        titlesForTable{end+1} = strcat(mediaTypesSelected{i},'vs',mediaTypesSelected{j});
        mediaIndsDict(i,j) = k;
        k= k+1;
    end
end

cancerTypesVec = [cancerTypes{1}',cancerTypes{2}'];
cancerTypesVec = cellfun(@(x) strip(x),cancerTypesVec,'un',0);
[cancerTypesVecGrouped,cancersType] = findgroups(cancerTypesVec);

cancerTypesIds = [];
for i=1:max(cancerTypesVecGrouped)
    if(length(find(cancerTypesVecGrouped==i))>10)
        cancerTypesIds =[cancerTypesIds,i];
    end
end

cancerTypesIds = sort(cancerTypesIds);
cancersType = cancersType(cancerTypesIds);

[ genesEssen,idsInDepMap,~ ] = FindGenesSigmaThresh( depMapData, thresh,0 );
[isExists,locs] = ismember(genesToCheck,genesEssen);
idsInDepMapSelected = idsInDepMap(locs(isExists));

pVals = zeros(length(idsInDepMapSelected),length(cancerTypesIds))+NaN;
rVals = zeros(length(idsInDepMapSelected),length(cancerTypesIds))+NaN;

for k=1:length(idsInDepMapSelected)
    index = idsInDepMapSelected(k);
    col = 1;
    for i=1:length(mediaTypesSelected)
        data1 = depMapData.data(index,mediaInds{i});
        for j=i+1:length(mediaTypesSelected)
            data2 = depMapData.data(index,mediaInds{j});
            if(~(sum(~isnan(data1))==0 || sum(~isnan(data2))==0))
                
                concData =[data1,data2]';
                isNotNanData1 = ~isnan(data1);
                isNotNanData2 = ~isnan(data2);
                mediaVec = [zeros(1,sum(isNotNanData1)),zeros(1,sum(isNotNanData2))+1]';
                isNotNan = ~isnan(concData);
                cancerVecAfterNans = cancerTypesVecGrouped(isNotNan)';
                concData = concData(isNotNan);
                
                for t = 1:length(cancerTypesIds)
                    cancerId = cancerTypesIds(t);
                    cancersVec = cancerVecAfterNans;
                    cancersVec(find(cancersVec~=cancerId)) = 0;
                    cancersVec(find(cancersVec==cancerId)) = 1;
                    [rho,p] = partialcorr(concData,mediaVec,cancersVec,'type','Spearman');
                    pVals(k,t) = p;
                    rVals(k,t) = rho;
                end                
            end
        end
    end
end
res.genes = depMapData.genes(idsInDepMapSelected);
res.pVals = pVals;
res.cancers = cancersType;
res.rVals = rVals;
end
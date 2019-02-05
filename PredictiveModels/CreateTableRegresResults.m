function [ T ] = CreateTableRegresResults(res)
numFeaturesToShow = 10;

genes = res.genes;
p = res.pValues;
r = res.rValues;
bestImpFeatures = res.bestImpFeatures;
bestImpFeaturesRank = res.bestImpFeaturesRank;
pFDR = res.pValuesFDR;


featureImp = cell(length(bestImpFeaturesRank),1);
featuresStr = cell(length(bestImpFeaturesRank),1);
for i=1:length(bestImpFeaturesRank)
    data = bestImpFeaturesRank{i};
    currentBesrFeatures = bestImpFeatures{i};
    numFeat = min(length(data),numFeaturesToShow);
    featureImp{i} = num2str(round(data(1:numFeat)*100)/100);
    featuresStr{i} = strjoin(currentBesrFeatures(1:numFeat),', ');
end
T1 = cell2table([genes,featuresStr,featureImp],'VariableNames',{'Gene','BestFeatures','FeatureImportance'});
T2 = array2table([p,r,pFDR],'VariableNames',{'p','r','pFDRWithRandom'});

T = [T1,T2];
end


function pValuesRandom = RandomPermutationPredictiveModels(numPerm,varTypes,genesSelected,...
    takeAllGenes,geneExp,cnv,mut,screen,screenMetadata,resDist,equalWeights)
pValuesRandom = zeros(numPerm,1);
randInds = randi(length(genesSelected),numPerm,1);
parfor i=1:numPerm
    if(mod(i,1000)==0)
        disp(i);
    end
    id = randInds(i);
    res = CreateDataForPredictor(genesSelected{id}, geneExp,cnv,...
        mut,screen,screenMetadata,takeAllGenes,resDist);

    features = find(contains(res.varNames,varTypes));
    
    X = res.data(:,features);
    Y = res.Z;
    Y = Y(randperm(length(Y)));
    getImportance = 0;
    [r,p,~,~] = GetPredictor(X,Y,res.isCategorized(features),getImportance,equalWeights);
    pValuesRandom(i)  = p;
end
end

function [pVal,countersArr,countGenesFound] = CompareToRandom(resToCheck,resExpEsnAll,numReps,boolUseBonferoni)

[countGenesFound,~] = FindGenesHavingCorrelation(resToCheck,boolUseBonferoni);
countersArr = zeros(numReps,1);
count = 0;
parfor i=1:numReps
    resExpEsnRand = CorrelationsEsnExpRandom(resToCheck,resExpEsnAll);
    [countGenesFoundRand,~] = FindGenesHavingCorrelation(resExpEsnRand,boolUseBonferoni);
    if(countGenesFoundRand>=countGenesFound)
        count = count+1;
    end
    countersArr(i) =countGenesFoundRand;
end
pVal= (1+count)/numReps;
end

function resExpEsnRand = CorrelationsEsnExpRandom(resExpEsn,resExpEsnAll)
resExpEsnRand = struct;
resExpEsnRand.genesEsn = resExpEsnAll.genesEsn;
resExpEsnRand.genesExp = resExpEsnAll.genesExp;

resExpEsnRand.r = zeros(size(resExpEsnAll.r))+NaN;
resExpEsnRand.p = zeros(size(resExpEsnAll.r))+NaN;

for i=1:length(resExpEsnRand.genesEsn)
    locsRelatedGenes = find(~isnan(resExpEsn.p(i,:)));
    if(~isempty(locsRelatedGenes)) 
        locsNotNanInAll = find(~isnan(resExpEsnAll.r(i,:)));
        p = randperm(length(locsNotNanInAll));
        p = p(1:length(locsRelatedGenes));
        locs = locsNotNanInAll(p);
        resExpEsnRand.r(i,locs) = resExpEsnAll.r(i,locs);
        resExpEsnRand.p(i,locs) = resExpEsnAll.p(i,locs);
    end
end

resExpEsnRand.pFDRBonferoni = calcFDR(resExpEsnRand.p,1);
resExpEsnRand.pFDR = calcFDR(resExpEsnRand.p,0);
end
addPaths;
configFile;
global Config;
maxCurrency = Config.CURRENCY_MET_THRESH;
maxNeighbors = Config.MAX_NUM_OF_NEIGHBORS;

% Load screens - rnaiScreen and crisprScreen structs- with the fields: 
% data, genes, celllines

% Load copy-number-variation(ccleCNV struct), gene expression(ccleExp) structs
% with the fields: common (gene symbol), caseId (cell-line), data 

[ genesSelectedAll,~ ] = FindGenesSigmaThresh(rnaiScreen,2,0);
load(Config.METABOLIC_GENES_PATH);
genesSelected = metGenes(ismember(metGenes,genesSelectedAll));

%load distances from each gene by the metabolic network
load(Config.METABOLIC_NETWORK_DISTANCES)

%Calculate correlation between essntiality and expression of all genes
%RNAi
[resExpEsn] = CorrelationsBetweenEsnToExp(rnaiScreen,ccleExp,genesSelected);
resExpEsnSelf = CorrelationsEsnExpIso(resExpEsn,resDist,-1,-1,1,Inf);
resExpEsnIso = CorrelationsEsnExpIso(resExpEsn,resDist,0,0,0,Inf);
resExpEsnRelated = CorrelationsEsnExpIso(resExpEsn,resDist,1,maxCurrency,0,maxNeighbors);

%Calculate correlation between essntiality and expression of all genes
%CRISPR
[exists,locs] = ismember(genesSelected,crisprScreen.genes);
genesSelectedCRISPR = genesSelected(exists);
[resExpEsnCRISPR] = CorrelationsBetweenEsnToExp(crisprScreen,ccleExp,genesSelectedCRISPR);
resExpEsnSelfCRISPR = CorrelationsEsnExpIso(resExpEsnCRISPR,resDist,-1,-1,1,Inf);
resExpEsnIsoCRISPR = CorrelationsEsnExpIso(resExpEsnCRISPR,resDist,0,0,0,Inf);
resExpEsnRelatedCRISPR = CorrelationsEsnExpIso(resExpEsnCRISPR,resDist,1,maxCurrency,0,maxNeighbors);

%Calculate correlation between essntiality and CNV of all genes
%RNAi
[resCnvEsn] = CorrelationsBetweenEsnToExp(rnaiScreen,ccleCNV,genesSelected);
resCnvEsnSelf = CorrelationsEsnExpIso(resCnvEsn,resDist,-1,-1,1,Inf);
resCnvEsnIso = CorrelationsEsnExpIso(resCnvEsn,resDist,0,0,0,Inf);
resCnvEsnRelated = CorrelationsEsnExpIso(resCnvEsn,resDist,1,maxCurrency,0,maxNeighbors);

%Calculate correlation between essntiality and CNV of all genes
%CRISPR
[resCnvEsnCRISPR] = CorrelationsBetweenEsnToExp(crisprScreen,ccleCNV,genesSelectedCRISPR);
resCnvEsnSelfCRISPR = CorrelationsEsnExpIso(resCnvEsnCRISPR,resDist,-1,-1,1,Inf);
resCnvEsnIsoCRISPR = CorrelationsEsnExpIso(resCnvEsnCRISPR,resDist,0,0,0,Inf);
resCnvEsnRelatedCRISPR = CorrelationsEsnExpIso(resCnvEsnCRISPR,resDist,1,maxCurrency,0,maxNeighbors);

numReps = Config.NUM_REP_RAND_CORRS;
[pValRelatedVsAll,countersRelatedVsAll,countRelated] = CompareToRandom(resExpEsnRelated,resExpEsn,numReps,0);
[pValSelfVsAll,countersSelfVsAll,countSelf] = CompareToRandom(resExpEsnSelf,resExpEsn,numReps,0);
[pValIsoVsAll,countersIsoVsAll,countISo] = CompareToRandom(resExpEsnIso,resExpEsnAll,numReps,0);
countsFound = [countSelf;countISo;countRelated]./length(genesSelected);
meansRand = [mean(countersSelfVsAll);mean(countersIsoVsAll);mean(countersRelatedVsAll)]./length(genesSelected);
stdRand = [std(countersSelfVsAll);std(countersIsoVsAll);std(countersRelatedVsAll)]./length(genesSelected);
PlotPairBarGraphWithStd(countsFound,meansRand,stdRand,{'Self','Isozymes','Metabolic Neighbors'},'2a')

[pValRelatedVsAllCRISPR,countersRelatedVsAll,countRelated] = CompareToRandom(resExpEsnRelatedCRISPR,resExpEsnCRISPR,numReps,0);
[pValSelfVsAllCRISPR,countersSelfVsAll,countSelf] = CompareToRandom(resExpEsnSelfCRISPR,resExpEsnCRISPR,numReps,0);
[pValIsoVsAllCRISPR,countersIsoVsAll,countISo] = CompareToRandom(resExpEsnIsoCRISPR,resExpEsnCRISPR,numReps,0);
countsFound = [countSelf;countISo;countRelated]./length(genesSelectedCRISPR);
meansRand = [mean(countersSelfVsAll);mean(countersIsoVsAll);mean(countersRelatedVsAll)]./length(genesSelectedCRISPR);
stdRand = [std(countersSelfVsAll);std(countersIsoVsAll);std(countersRelatedVsAll)]./length(genesSelectedCRISPR);
PlotPairBarGraphWithStd(countsFound,meansRand,stdRand,{'Self','Isozymes','Metabolic Neighbors'},'2b')


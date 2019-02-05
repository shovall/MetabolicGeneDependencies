addPaths;
configFile;
global Config;
numRep = Config.NUM_REP_RANDOM_P_VALUES;

% Load screens - rnaiScreen and crisprScreen structs- with the fields: 
% data, genes, celllines

% Load screens metadata - rnaiMetadata and crisprMetadata - with the fields:
% cell lines, mediaType, primaryDisease, mediaTypeGrouped

% Load copy-number-variation(ccleCNV struct), gene expression(ccleExp) and
% oncogene/tumor suppressor(driverMut) struct with the fields: common (gene
% symbol), caseId (cell-line), data 

[ genesSelectedAll,~ ] = FindGenesSigmaThresh(rnaiScreen,2,0);
load(Config.METABOLIC_GENES_PATH);
genesSelected = metGenes(ismember(metGenes,genesSelectedAll));

driverMut = AddMutDataDiscrete(driverMut);
%load distances from each gene by the metabolic network
load(Config.METABOLIC_NETWORK_DISTANCES)

%%%RNAi
% Random permutation p-values
pValuesRandomNoMediaRNAi = RandomPermutationPredictiveModels(numRep,{'Current','Related','Mut'},genesSelected,...
    0,ccleExp,ccleCNV,driverMut,rnaiScreen,rnaiMetadata,resDist,0);

pValuesRandomRNAi = RandomPermutationPredictiveModels(numRep,{'Current','Related','Mut','Media'},genesSelected,...
    0,ccleExp,ccleCNV,driverMut,rnaiScreen,rnaiMetadata,resDist,0);

pValuesRandomRNAiCancer = RandomPermutationPredictiveModels(numRep,{'Current','Related','Mut','Media','CancerType'},genesSelected,...
    0,ccleExp,ccleCNV,driverMut,rnaiScreen,rnaiMetadata,resDist,0);

% Predictors calculation
[resRelatedRNAiNoMedia] = GetPredictorsForGenes(genesSelected ,ccleExp,ccleCNV, driverMut,rnaiScreen,rnaiMetadata,{'Current','Related','Mut'},0,resDist,0);
[resRelatedRNAiNoMedia.pValuesFDR,resRelatedRNAiNoMedia.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(resRelatedRNAiNoMedia.pValues,pValuesRandomNoMediaRNAi);

[resRelatedRNAi] = GetPredictorsForGenes(genesSelected ,ccleExp,ccleCNV, driverMut,rnaiScreen,rnaiMetadata,{'Current','Related','Mut','Media'},0,resDist,0);
[resRelatedRNAi.pValuesFDR,resRelatedRNAi.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(resRelatedRNAi.pValues,pValuesRandomRNAi);

[resRelatedRNAiCancer] = GetPredictorsForGenes(genesSelected ,ccleExp,ccleCNV, driverMut,rnaiScreen,rnaiMetadata,{'Current','Related','Mut','Media','CancerType'},0,resDist,0);
[resRelatedRNAiCancer.pValuesFDR,resRelatedRNAiCancer.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(resRelatedRNAiCancer.pValues,pValuesRandomRNAiCancer);

%%% CRISPR---
% Random permutation p-values
[exists,locs] = ismember(genesSelected,crisprScreen.genes);
genesSelectedCERES = genesSelected(exists);
pValuesRandomNoMediaCRISPR = RandomPermutationPredictiveModels(numRep,{'Current','Related','Mut'},genesSelectedCERES,...
    0,ccleExp,ccleCNV,driverMut,crisprScreen,crisprMetadata,resDist,1);

pValuesRandomCRISPR = RandomPermutationPredictiveModels(numRep,{'Current','Related','Mut','Media'},genesSelectedCERES,...
    0,ccleExp,ccleCNV,driverMut,crisprScreen,crisprMetadata,resDist,1);

pValuesRandomCRISPRCancer = RandomPermutationPredictiveModels(numRep,{'Current','Related','Mut','Media','CancerType'},genesSelectedCERES,...
    0,ccleExp,ccleCNV,driverMut,crisprScreen,crisprMetadata,resDist,1);

% Predictors calculation
[resRelatedCRISPRNoMedia] = GetPredictorsForGenes(genesSelectedCERES ,ccleExp,ccleCNV, driverMut,crisprScreen,crisprMetadata,{'Current','Related','Mut'},0,resDist,1);
[resRelatedCRISPRNoMedia.pValuesFDR,resRelatedCRISPRNoMedia.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(resRelatedCRISPRNoMedia.pValues,pValuesRandomNoMediaCRISPR);

[resRelatedCRISPR] = GetPredictorsForGenes(genesSelectedCERES ,ccleExp,ccleCNV, driverMut,crisprScreen,crisprMetadata,{'Current','Related','Mut','Media'},0,resDist,1);
[resRelatedCRISPR.pValuesFDR,resRelatedCRISPR.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(resRelatedCRISPR.pValues,pValuesRandomCRISPR);

[resRelatedCRISPRCancer] = GetPredictorsForGenes(genesSelectedCERES ,ccleExp,ccleCNV, driverMut,crisprScreen,crisprMetadata,{'Current','Related','Mut','Media','CancerType'},0,resDist,1);
[resRelatedCRISPRCancer.pValuesFDR,resRelatedCRISPRCancer.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(resRelatedCRISPRCancer.pValues,pValuesRandomCRISPRCancer);

%RNAi - Compare models to random neighbors 
reps = Config.NUM_REP_RAND_NEIGHBORS;
randNeighborsModelsRNAi = cell(reps,1);
fracsRNAi = zeros(reps,1);
for i=1:reps
    cur = GetPredictorsForGenes(genesSelected ,ccleExp,ccleCNV, driverMut,rnaiScreen,rnaiMetadata,{'Current','Related','Mut'},0,resDist,0,1);
    [cur.pValuesFDR,cur.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(cur.pValues,pValuesRandomNoMediaRNAi);
    randNeighborsModelsRNAi{i} = cur;
    fracsRNAi(i) = length(find(randNeighborsModelsRNAi{i}.pValuesFDR<0.05))/length(randNeighborsModelsRNAi{i}.pValuesFDR);
end

%CRISPR - Compare models to random neighbors 
randNeighborsModelsCRISPR = cell(reps,1);
fracs = zeros(reps,1);
for i=1:reps
    cur = GetPredictorsForGenes(genesSelectedCERES,ccleExp,ccleCNV, driverMut,crisprScreen,crisprMetadata,{'Current','Related','Mut'},0,resDist,1,1);
    [cur.pValuesFDR,cur.pValuesWithRand] = CalcFDRWithOtherRandPermPVals(cur.pValues,pValuesRandomNoMediaCRISPR);
    randNeighborsModelsCRISPR{i} = cur;
    fracs(i) = length(find(randNeighborsModelsCRISPR{i}.pValuesFDR<0.05))/length(randNeighborsModelsCRISPR{i}.pValuesFDR);
end

[ resRelatedRNAiCancer.T ] = CreateTableRegresResults(resRelatedRNAiCancer);
writetable(resRelatedRNAiCancer.T,'RNAi models.xlsx');

[ resRelatedCRISPRCancer.T ] = CreateTableRegresResults(resRelatedCRISPRCancer);
writetable(resRelatedCRISPRCancer.T,'CRISPR models.xlsx');

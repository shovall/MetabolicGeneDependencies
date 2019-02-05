addPaths;
configFile;
global Config;
screenType = 'RNAi'; %Or CRISPR
analysisTitle = 'RNAi6'; % Will be used as the prefix for the histograms plotted
thresh = 6; % Select 6-sigma genes

% Load screens - rnaiScreen and crisprScreen structs- with the fields: 
% data, genes, celllines
% Load screens metadata - rnaiMetadata and crisprMetadata - with the fields:
% cell lines, mediaType, primaryDisease, mediaTypeGrouped, conditions
load(Config.METABOLIC_GENES_PATH);

if(screenType=="RNAi")
    screen = rnaiScreen;
    screenMetadata = rnaiMetadata;
else
    screen = rnaiScreen;
    screenMetadata = rnaiMetdata;
end

load(Config.METABOLIC_GENES_PATH);
[metGenes,selectedMetGenes,selectedNonMet] = GetSelectedGenes(rnaiScreen,metGenes,thresh);

% Selct cell lines cultured in DMEM/RPMI
DMEMCells = screenMetadata.celllines(contains(screenMetadata.mediaType,'DMEM'));
RPMICells = screenMetadata.celllines(contains(screenMetadata.mediaType,'RPMI'));
mediaTypesCells = {DMEMCells,RPMICells};
mediaTypesSelected = {'DMEM','RPMI'};

plotHistograms = false;
% Compare distributions of dependency scores between cell lines cultured in
% different media types
%For metabolic genes
[genesRes,pVals,pValsFDR,pSpearPartCond,r,T] = GetTableMediaScreening(plotHistograms,analysisTitle,selectedMetGenes,screen,screenMetadata,mediaTypesCells,mediaTypesSelected,0);
% For non-metabolic genes
[genesResAll,pValsAll,pValsFDRAll,pSpearPartCondAll,rAll,TAll] = GetTableMediaScreening(plotHistograms,strcat(analysisTitle,'all'),selectedNonMet,screen,screenMetadata,mediaTypesCells,mediaTypesSelected,0);

% FDR correction for both metabolic and non-metabolic genes analyzed
pMetAndAll = [pVals;pValsAll];
pMetAndAllFDR = CalcFDR(pMetAndAll);

%Create figures
createHistForPValues;
CreateBarGraphMostSensitiveGenes(pVals,pMetAndAllFDR,genesRes,r,'2e');

% Controlling for cancer type (using Sperman partial correlation)
[resCancerCancel] = GetTableMediaScreeningCancerCancel(selectedMetGenes,screen,screenMetadata,mediaTypesCells,mediaTypesSelected,0);

function [metGenes,selectedMetGenes,selectedNonMet] = GetSelectedGenes(depMapData,metGenes,thresh)
[ genesThresh,~,~ ] = FindGenesSigmaThresh(depMapData,thresh,0);
selectedMetGenes = metGenes(ismember(metGenes,genesThresh));
[exists,~] =ismember(depMapData.genes,metGenes);
selectedNonMet = depMapData.genes(~exists);
selectedNonMet = intersect(selectedNonMet,genesThresh);
end
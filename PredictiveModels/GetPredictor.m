function [r,p,featSelected,imp,Ensemble] = GetPredictor(X,Y,isCategorized,getImportance,equalWeights)
catVec = find(isCategorized(:));
if(nargin>4 && equalWeights==1)
    weights = ones(size(Y));
else
    weights = GetWeigthsForPredictor(Y);
end
if(getImportance)
    impStr = 'On';
else
    impStr = 'Off';
end

if(isempty(X))
    r = NaN;
    p= NaN;
    featSelected = NaN;
    imp= NaN;
    Ensemble = NaN;
    return;
end

Ensemble = GetEnsemble(X,Y,catVec,weights,impStr);
YHat = Ensemble.oobPredict;
[r,p] = corr(YHat,Y);

if (~getImportance)
    featSelected = NaN;
    imp = NaN;
    return;
end
impAll = Ensemble.OOBPermutedPredictorDeltaError;
featuresImportant = find(impAll>0);
featSelected = featuresImportant;
catVec = find(isCategorized(featSelected));
if(isempty(featuresImportant))
    imp = [];
else
    Ensemble2 = GetEnsemble(X(:,featuresImportant),Y,catVec,weights,impStr);
    imp = Ensemble2.OOBPermutedPredictorDeltaError;
end
end

function Ensemble = GetEnsemble(X,Y,catVec,weights,impStr)
    global Config;
    Ensemble = TreeBagger(Config.NUM_TREES,X,Y,'MinLeafSize',Config.MIN_LEAF_SIZE,...
    'Method','regression',...
    'CategoricalPredictors',catVec,...
    'Weights',weights,'OOBPrediction','On',...
    'OOBPredictorImportance',impStr);
end
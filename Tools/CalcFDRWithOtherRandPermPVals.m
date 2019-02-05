function [pValuesFDR,pValuesWithRand] = CalcFDRWithOtherRandPermPVals(pValues,pValuesRandPerm)

pValuesWithRand = zeros(size(pValues));

for i=1:length(pValuesWithRand)
   pValuesWithRand(i) = (1+length(find(pValuesRandPerm<pValues(i))))/(1+length(find(pValuesRandPerm>=0)));
end

pValuesFDR = CalcFDR(pValuesWithRand);
end
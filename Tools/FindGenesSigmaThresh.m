function [ genes,ids,amounts ] = FindGenesSigmaThresh( depMapData,thresh, threshNumCellLineEssen)

ids = [];
amounts = [];

for i=1:size(depMapData.data)
    amount =length(find(depMapData.data(i,:)<=-thresh));
    if(amount>threshNumCellLineEssen)
        ids(end+1) = i;
        amounts(end+1) = amount;
    end
end

genes = depMapData.genes(ids);
end


function weights = GetWeigthsForPredictor(Y)
maxVal= 0.05;
maxVal = maxVal*length(Y);
idsEssen = find(Y<-2);
weights = zeros(size(Y))+1;
weights(idsEssen) = (length(weights)-length(idsEssen))/length(idsEssen);

weights = min(weights,maxVal);
weights = weights./sum(weights);
end

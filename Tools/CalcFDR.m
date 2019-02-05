function [ adjustedP ] = CalcFDR( p, doBonferoni )
if(nargin<2)
    doBonferoni=0;
end

origSize = size(p);
p_as_1D_vector = reshape(p,[(size(p,2)*size(p,1)),1]);

if(~doBonferoni)
    adjustedFDR = mafdr(p_as_1D_vector,'BHFDR',true);
else
    numElements = sum(~isnan(p_as_1D_vector));
    adjustedFDR = p_as_1D_vector*numElements;
    adjustedFDR(~isnan(adjustedFDR)) = min(1,adjustedFDR(~isnan(adjustedFDR)));
end
adjustedP = reshape(adjustedFDR,origSize);

end


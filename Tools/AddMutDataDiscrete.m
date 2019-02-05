function [ mutDataStruct ] = AddMutDataDiscrete( mutDataStruct )

discrete = zeros(size(mutDataStruct.data));

for i=1: size(mutDataStruct.data,1)
    for j=1:size(mutDataStruct.data,2)
        if strcmp(mutDataStruct.data(i,j),'NaN')
            discrete(i,j) = 0;
        else
            discrete(i,j) = 1;
        end
    end 
end

mutDataStruct.discreteData = discrete;

end


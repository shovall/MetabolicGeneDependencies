folders = {'Data','MediaEffect','MetabolicNeighbors','PredictiveModels','Tools'};

for i=1:length(folders)
    addpath(genpath(folders{i}));
end
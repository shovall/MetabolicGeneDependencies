global Config 
Config = {};

% Paths Only:
Config.FIGURES_FOLDER = 'C:\MediaScreeningProject\Figures';
Config.MEDIA_SCREENING_HISTOGRAMS = 'C:\MediaScreeningProject\Histograms'; 
Config.METABOLIC_GENES_PATH = 'Data/MetabolicGenes';
Config.METABOLIC_NETWORK_DISTANCES = 'Data/resDist';

fields = fieldnames(Config);
for i=1:numel(fields)
  Config.(fields{i}) = GetFullPath(Config.(fields{i}));
end

% Not Paths
Config.NUM_REP_RAND_CORRS = 1000;
Config.NUM_REP_RANDOM_P_VALUES = 25000;
Config.NUM_REP_RAND_NEIGHBORS = 100;
Config.NUM_TREES = 100;
Config.MIN_LEAF_SIZE = 5;
Config.CANCER_LINEAGE_MIN_NUMBER = 10; 
Config.MAX_NUM_OF_NEIGHBORS = 50;
Config.CURRENCY_MET_THRESH = 50;
Config.FRAC_NAN_IGNORE_FEATURE = 0.5;

clear('fields','i');
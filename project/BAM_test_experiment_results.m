

% BAM_test_experiment_results

% needed variables...
% ====================
% - anTestData - (dim x num_test_points)
% - anTrainingData - (dim x num_training_points)
% - anTestGenreIndex  - (1 x num_test_points,  values 1-6)
% - anTrainingGenreIndex  - (1 x num_training_points, values 1-6)

% IF YOU NEED TO TEST...
% The file below in the 'load' command contains the variables above, 
% which is all you need.  Then hit 'run' and see results.
%
%>> load(strcat(sDir,'BAM_test_experiment_results_testdata.mat'));




% setting knn search parameters
nNumNbrs = 7; 
test_equal_training_data = 1;   % set to 1 if anTestData is same as anTrainingData
                                % otherwise set to 0

% generating the knn data info
[anTestData_nbr_index anTestData_nbr_dist] = ...
    BAM_get_knn_data_info(anTrainingData, anTestData, nNumNbrs, test_equal_training_data);


% performing any adjustments on the 'raw' nearest neighbor info...
% (e.g. removing neighbors known to be poor, or removing all of a given
%  genre if there are below some threshold, etc.)



% setting the knn distance distribution parameters
% (see 'BAM_make_knn_data_dist()' code for how they affect neighbor vote
% weighting)
nNumLabel = 6;  % number of labels (genres)
DistanceWeightOn = 1;
DistanceWeightType = 1;
DistanceWeightFactor = 1;
DistanceFilterThreshold = Inf;
LabelWeightOn = 1;  
LabelWeightType = 2; 
LabelWeightFactor = 0.075;


% making the knn data distribution and genre probability
[anTestData_nbr_dist_wt anTestData_label_prob] = BAM_make_knn_data_dist(...
    anTestData_nbr_index, anTestData_nbr_dist, nNumLabel, anTrainingGenreIndex, ...
    DistanceWeightOn, DistanceWeightType, DistanceWeightFactor, DistanceFilterThreshold, ...
    LabelWeightOn, LabelWeightType, LabelWeightFactor);


% performing any further adjustments on the knn data distribution and genre
% probability



% setting the default label index
DefaultLabelIndex = 6;



% performing the label (genre) assignment based on knn data distribution
% and genre probability
[anTestGenreIndex_assigned] = BAM_assign_knn_dist_label(anTestData_label_prob, ...
        anTestData_nbr_index, anTestData_nbr_dist, anTestData_nbr_dist_wt, DefaultLabelIndex);
    

% producing the performance results - confusion matrix, percentage correct    
[anConfusionMatrix, anPercentCorrect, nAllPercentCorrect] = ...
    BAM_make_confusion_matrix(anTestGenreIndex_assigned, anTestGenreIndex, nNumLabel)




% IDEA!!! - when performing cross-validation or unknown test data,
%           first perform algorithm on TrainingData points,
%           then somehow weight 'correctly classified songs'
%           preferentially to 'incorrectly classified songs'????




%{
Scripts to do the following:
+ Read the emebeddedDataPoints as Row-Vectors
+ Determine number of iterations over which Clustering should be
implemented.
+ Partition into TRAIN-TEST datasets. 
+ Compute the confusionMatrix for each iteration.
+ Determine the Mean and StdDev of ConfusionMatrices.

Plots:
1. Imagesc of the Mean and Covariance ConfusionMatrix for the Test-Data
(Voting of the ClusterMemebership of TestData by the TrainingData.)

%}


% ========================================================================
% Determine the Training and Test Data: 

function [mean_ConfusionMatrix, stdDev_ConfusionMatrix, mean_GenrePercentCorrect, std_GenrePercentCorrect, avgSuccessRate, stdDev_SuccessRate] =  run_experiment_train_test_BAM(tracksDirName, dataSet, totalIter)
%{
Input: 
tracksDirName:  Directory containing tracks
dataSet:        Coordinates of songs as row-vectors, in the Embedded Space.
totalIter:      Number of iterations of random training and test/data to be
created.

Output:
mean_ConfusionMatrix:      mean ConfusionMatrix for the Genres
stdDev_ConfusionMatrix:     stddev of ConfusionMatrix for the Genres

%}

% Parameters for the code to run.
% totalIter = 30;
% tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';

% -------------------------------------------------------------------
% Iterate the following section, to compute the ConfusionMatrix for given
% Pipeline. 

[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

allConfusionMatrices = zeros(length(genreKeys),length(genreKeys),totalIter);
allGenrePercentCorrect = zeros(totalIter,length(genreKeys));
totalPercentCorrect = zeros(totalIter,1);
    
for iter = 1:totalIter
    iter
    % -------------------------------------------------------------------
    % Create a dictionary of indices for each genreType.
    
    trainData_Indices = [];
    trainData_genreID = [];
    testData_Indices = [];
    testData_genreID = [];

    
    % for each genre, identify 80% as Training, 20% as Test.
    for genreIndx = 1:length(genreKeys)
        songsOfAGenre = find(songGenres == genreIndx); % indices of all songs belonging to
        
        % randomly permute the index-vector
        totalSongsInGenre = length(songsOfAGenre);
        permutation = randperm(totalSongsInGenre);
        songsOfAGenre = songsOfAGenre(permutation);
        
        % extract indices of the Train and Test Data
        lengthTrainData = ceil(0.8*totalSongsInGenre);
        
        trainData_Indx = songsOfAGenre(1:lengthTrainData);
        trainData_genre = genreIndx*ones(1,length(trainData_Indx));
        
        testData_Indx = songsOfAGenre(lengthTrainData+1: end);
        testData_genre = genreIndx*ones(1,length(testData_Indx));
        
        % store/archive into the main vector
        trainData_Indices = [trainData_Indices trainData_Indx];
        trainData_genreID = [trainData_genreID trainData_genre];
        
        testData_Indices = [testData_Indices testData_Indx];
        testData_genreID = [testData_genreID testData_genre];
        
        
        %-------------------------------------------------------------
        % Update the Dictionary:
        % insert into container! - for genreSpecific archival
        if genreIndx ~= 1
            tmpTrainData = containers.Map(genreIndx, [trainData_Indx]);
            tmpTestData = containers.Map(genreIndx, [testData_Indx]);
            
            trainData_dict = [trainData_dict; tmpTrainData];
            testData_dict = [testData_dict; tmpTestData];
            %         trainData(genreIndx) = trainData_Indx;
            %         testData(genreIndx) = testData_Indx;
        else
            trainData_dict = containers.Map(genreIndx, [trainData_Indx]);
            testData_dict = containers.Map(genreIndx, [testData_Indx]);
        end
    end
    
    %     trainData_dict.keys()
    %     trainData_dict.values()
    
    % Determine Size of Training and Testing Datasets!
    %     size(trainData_Indices)
    %     size(testData_Indices)
    % =========================================================================
    % Extract the Training and Testing Datasets from the givenDataset
    
    embeddedPoints_train = dataSet(trainData_Indices,:);
    trainData_genreID;
    
    embeddedPoints_test = dataSet(testData_Indices,:);
    testData_genreID;
    % ========================================================================
    % BAM_Cluster_Membership
    
    
% input variables...
%  - anTestData (d x num_test_samples)
%  - anTrainingData (d x num_training_samples)
%    (both above are the data vectors in whatever dimensions desired)
%
%  - anTestGenreIndex (1 x num_test_samples)
%  - anTrainingGenreIndex (1 x num_training_samples)
%    (both above are lists of genre index (1-6)
%
%  PARAMETERS...
%  - num_nbrs = the number of neighbors
%  - weight_type = (1-3) for distance weighting type
%  - weight_sigma = for Gaussian distance weighting
%  - test_equal_training_data = set to 1 if test data is the same as
%                               training data
%
% OUTPUTS...
%  - anTestData_ConfusionMatrix
%  - anTestData_PercentCorrect
%  - anTestData_AllPercentCorrect


anTrainingData = embeddedPoints_train';
anTestData = embeddedPoints_test';

anTestGenreIndex = reshape(testData_genreID,1,[]);
anTrainingGenreIndex = reshape(trainData_genreID,1,[]);

BAM_test_experiment_results

% anPercentCorrect, nAllPercentCorrect

% confusionMatrix = anTestData_ConfusionMatrix(:,:,4);
sumConf = sum(anConfusionMatrix,1);
factorMatrix = repmat(sumConf,length(anConfusionMatrix),1);
confusionMatrix = anConfusionMatrix./factorMatrix;

% ========================================================================

    allConfusionMatrices(:,:,iter) = confusionMatrix;
    allGenrePercentCorrect(1,:,iter) = anPercentCorrect;
    totalPercentCorrect(iter,1) = nAllPercentCorrect;
end

% Computing Mean and Std-Dev of the Confusion Matrices
mean_ConfusionMatrix = 1/totalIter * sum(allConfusionMatrices,3)

% Computing the STD-DEV ConfusionMatrix
mean_ConfusionMatrix_REPMAT = repmat(mean_ConfusionMatrix,1,1,totalIter);
size(mean_ConfusionMatrix_REPMAT)
size(allConfusionMatrices)
var_ConfusionMatrix = 1/totalIter * sum((allConfusionMatrices - mean_ConfusionMatrix_REPMAT).^2,3);

stdDev_ConfusionMatrix = var_ConfusionMatrix.^(0.5);

% computing the totalAverage and stdDev
mean_GenrePercentCorrect = diag(mean_ConfusionMatrix); 
std_GenrePercentCorrect = diag(stdDev_ConfusionMatrix);

avgSuccessRate = sum(totalPercentCorrect)/totalIter;
mu_vec = repmat(avgSuccessRate,totalIter,1);
stdDev_SuccessRate = sqrt(sum((totalPercentCorrect - mu_vec).^2)/totalIter);


end


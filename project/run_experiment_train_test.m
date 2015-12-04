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

function [mean_ConfusionMatrix, stdDev_ConfusionMatrix] =  run_experiment_train_test(tracksDirName, dataSet, totalIter)
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
    % Determine Cluster Memebership for each of the test-songs
    
    % clusterMembership_Test
    testData_genreID_pred = [];
    for testsample = 1:length(testData_Indices)
        % pointwise determine the cluterMembership for each testSong.
        % Reads the PROB-VECTOR of ClusetrMembership, takes the
        % MaximumValue-Membership as its GenreID_Predicted.
        sampleCoord = embeddedPoints_test(testsample,:);
        sampleCoord_genreID_true = testData_genreID(testsample);
        clusterMembershipProb = getClusterMembership(sampleCoord,embeddedPoints_train,trainData_genreID);
        [M, sampleCoord_genreID_pred] = max(clusterMembershipProb);
        testData_genreID_pred(testsample) = sampleCoord_genreID_pred;
    end
    
    % determine total number of Genres
    totalGenres = length(genreKeys);
    
    % intialize ConfusionMatrix for given iteration
    confusionMatrix = zeros(totalGenres,totalGenres);
    
    for genreID = 1:totalGenres
        % identify indices of songs belonging to a genre, 
        % within TRUTH-DATSET of the Training Dataset.
        SongsIndxOfGivenGenre = (testData_genreID == genreID);
        totalSongsOfGenre = sum(SongsIndxOfGivenGenre);
        
        % Determine PREDICTION-ACCURACY for songs WITHIN GENRES.
        % search within given subset of test_songs:
        Pred_SongsIndxOfGivenGenre = testData_genreID_pred.*SongsIndxOfGivenGenre;
        
        for j=1:totalGenres
            % total Songs estimated to belong to particular genre
            numSongsEstAsGivenGenre = sum(Pred_SongsIndxOfGivenGenre == j);
            confusionMatrix(genreID,j) = numSongsEstAsGivenGenre/totalSongsOfGenre;
        end
    end
    
    allConfusionMatrices(:,:,iter) = confusionMatrix;
end

% Computing Mean and Std-Dev of the Confusion Matrices
mean_ConfusionMatrix = 1/totalIter * sum(allConfusionMatrices,3)

% Computing the STD-DEV ConfusionMatrix
mean_ConfusionMatrix_REPMAT = repmat(mean_ConfusionMatrix,1,1,totalIter);
var_ConfusionMatrix = 1/totalIter * sum((allConfusionMatrices - mean_ConfusionMatrix_REPMAT).^2,3);

stdDev_ConfusionMatrix = var_ConfusionMatrix.^(0.5);

%{
% Compute STD_DEV iteratively.
diffSum = zeros(size(mean_ConfusionMatrix));
for i = 1: totalIter
    diffSum = diffSum + (allConfusionMatrices(:,:,i) - mean_ConfusionMatrix).^2;
end
stdDev_ConfusionMatrix = (1/totalIter*diffSum).^(1/2)
%}


fig1 = figure(1)
colormap parula
subplot(2,1,1)
imagesc(mean_ConfusionMatrix);
title('MEAN Confusion Matrix')
colorbar

subplot(2,1,2)
imagesc(stdDev_ConfusionMatrix);
title('STD-DEV Confusion Matrix')
colorbar

end


%{

% ========================================================================
    % -------------------------------------------------------------------
    % Extract Distance Matrix for the Training and Testing Data
    
    mutualDistance_train = mutualDistance(trainData_Indices,trainData_Indices);
    mutualDistance_test = mutualDistance(testData_Indices,testData_Indices);
    
    % ========================================================================
    % refinedGraphEmbedding
    dimReducedSpace = 3;    % final dimension of reduced space
    embeddedPoints_train = refinedGraphEmbedding_distanceMatrix(mutualDistance_train, dimReducedSpace);
    
    embeddedPoints_test = refinedGraphEmbedding_distanceMatrix(mutualDistance_test, dimReducedSpace);
    % ========================================================================
    
    fig1 = figure(1)
    g_X = embeddedPoints_train(:,1);
    g_Y = embeddedPoints_train(:,2);
    g_Z = embeddedPoints_train(:,3);
    scatter3(g_X, g_Y, g_Z,'g')
    hold on
    g_X = embeddedPoints_test(:,1);
    g_Y = embeddedPoints_test(:,2);
    g_Z = embeddedPoints_test(:,3);
    scatter3(g_X, g_Y, g_Z,'r')
    title('Training and Test Data Projected')
    
    % -------------------------------------------------------------------
    % Generate Colored Plots of songs from each category
    
    close all
    fig2 = figure(2)
    % colormap(jet(length(genreKeys)));
    colormap(parula(length(genreKeys)));
    subplot(2,1,1)
    for i = 1:length(genreKeys)
        % determine indices of songs in training data belonging to  particular
        % genre.
        songsOfAGenre = find(trainData_genreID == i);
        hue = (length(genreKeys)+1-i)*ones(length(songsOfAGenre),1);
        rad = 30*ones(length(songsOfAGenre),1);
        g_X = embeddedPoints_train(songsOfAGenre,1);
        g_Y = embeddedPoints_train(songsOfAGenre,2);
        g_Z = embeddedPoints_train(songsOfAGenre,3);
        h = scatter3(g_X, g_Y, g_Z,rad,hue,'filled')
        hold on
    end
    legend(genreKeys)
    hold off
    
    subplot(2,1,2)
    for i = 1:length(genreKeys)
        % determine indices of songs in training data belonging to  particular
        % genre.
        songsOfAGenre = find(testData_genreID == i);
        hue = (length(genreKeys)+1-i)*ones(length(songsOfAGenre),1);
        rad = 30*ones(length(songsOfAGenre),1);
        g_X = embeddedPoints_test(songsOfAGenre,1);
        g_Y = embeddedPoints_test(songsOfAGenre,2);
        g_Z = embeddedPoints_test(songsOfAGenre,3);
        h = scatter3(g_X, g_Y, g_Z,rad,hue,'filled')
        hold on
    end
    legend(genreKeys)
    hold off
    % ---------------------------------------------------

% =========================================================================

%}

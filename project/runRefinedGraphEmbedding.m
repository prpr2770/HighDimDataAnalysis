%{

27 Novemeber!

Generate Graph Embedding of the entire database. 
Perform clustering of only the training Data. 
Then detemine the belongingness of the test Data. 

% ---------------------------------------------------------------
Script that does the following:
1. Read DistanceMatrix from folder. 
2. Perform RefinedGraph Embedding
3. Plot 3D vectors.


Part 2:
1. Determine clusters between the points
2. Identify the indexes of points in each cluster
3. Determine true Genre Classfication of songs, and compare with the
Clusters they belong to. 
4. What sort of plots are you doing?


%}

clear all;close all; clc
% dir containing tracks
tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\'
% dir with song features
g1c_SongFeatures_Dir = fullfile(tracksDirName ,'g1c_SongFeatures\')

% distance matrix files
distMatrix_txt = 'distance_matrix.txt';
distMatrix_mat = 'distance_matrix.mat';

% -------------------------------------------------------------------
% reading the distMatrix.mat file
matfile_distMatrix = fullfile(g1c_SongFeatures_Dir, distMatrix_mat);
distance_matrix = load(matfile_distMatrix);

mutualDistance = distance_matrix.D;

% -------------------------------------------------------------------
% refinedGraphEmbedding
dimReducedSpace = 3;    % final dimension of reduced space
embeddedPoints = refinedGraphEmbedding_distanceMatrix(mutualDistance, dimReducedSpace);

fig1 = figure(1)
g_X = embeddedPoints(:,1);
g_Y = embeddedPoints(:,2);
g_Z = embeddedPoints(:,3);
scatter3(g_X, g_Y, g_Z)


% -------------------------------------------------------------------
% Generate Colored Plots of songs from each category

% read the genreFile to determine the index 
songGenres_filename = 'songGenres.mat';
songGenres_full_filename = fullfile(tracksDirName,'g1c_SongFeatures\', songGenres_filename) ;
% save(songGenres_full_filename,'songGenres','genreKeys','wavfileGenreDict');
songGenresData = load(songGenres_full_filename);

songGenres = songGenresData.songGenres;
genreKeys = songGenresData.genreKeys;

totalColors = length(genreKeys);
colors = rand(totalColors,3);


% close all
fig2 = figure(2)
% colormap(jet(length(genreKeys)));
colormap(parula(length(genreKeys)));

for i = 1:length(genreKeys)
    songsOfAGenre = find(songGenres == i);
    hue = (length(genreKeys)+1-i)*ones(length(songsOfAGenre),1);
    rad = 30*ones(length(songsOfAGenre),1);
    g_X = embeddedPoints(songsOfAGenre,1);
    g_Y = embeddedPoints(songsOfAGenre,2);
    g_Z = embeddedPoints(songsOfAGenre,3);
    h = scatter3(g_X, g_Y, g_Z,rad,hue,'filled')
    hold on
end
legend(genreKeys)
hold off

% ========================================================================
% Clustering and Confusion Matrix Generation

% We have all the songs of database projected into lower dimensions. 
% We need to create the 80-20 partition of the dataset. 
% And then determine the membership of the songs in Training Data. 
% Further, compute the Confusion-Matrix. 

totalIter = 20;
allConfusionMatrices = zeros(length(genreKeys),length(genreKeys),totalIter);

for iter = 1:totalIter

    
    songGenres % Contains enumeration of the genreIDs for each song. 
    
    % -------------------------------------------------------------------
    % Create a dictionary of indices for each genreType.
    
    trainData_Indices = [];
    trainData_genreID = []; %vector stating genreID for each song represented by trainData_Indices.
    
    testData_Indices = [];
    testData_genreID = []; %vector stating genreID for each song represented by testData_Indices.
    
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
    
% %     trainData_dict.keys()
% %     trainData_dict.values()
% %     
% % %     Determine Size of Training and Testing Datasets!
% %     size(trainData_Indices)
% %     size(testData_Indices)
    % -------------------------------------------------------------------
    % Extract Coordinates of the training_data and test_data

    embeddedPoints_train = embeddedPoints(trainData_Indices,:);
    embeddedPoints_test = embeddedPoints(testData_Indices,:);
    
    % -------------------------------------------------------------------
    % clusterMembership_Test: Assigning a Cluster to each Test Song
    testData_genreID_est = []; % vector containing predicted-ClusterMembership for each Test-Song
    for testsample = 1:length(testData_Indices)
        sampleCoord = embeddedPoints_test(testsample,:);
        sampleCoord_genreID_true = testData_genreID(testsample);
        
        clusterMembershipProb = getClusterMembership(sampleCoord,embeddedPoints_train,trainData_genreID);
        
        % obtaining Single-Cluster Membership from the probability vector
        [M, sampleCoord_genreID_est] = max(clusterMembershipProb);
        testData_genreID_est(testsample) = sampleCoord_genreID_est;
    end
    
    
    % ---------------------------------------------------------------------
    totalGenres = length(genreKeys);
    confusionMatrix = zeros(totalGenres,totalGenres);
    
    for genreID = 1:totalGenres
        % identify indices of songs belonging to a genre
        SongsIndxOfGivenGenre = (testData_genreID == genreID);
        totalSongsOfGenre = sum(SongsIndxOfGivenGenre);
        
        %search within given subset of test_songs
        Est_SongsIndxOfGivenGenre = testData_genreID_est.*SongsIndxOfGivenGenre;
        
        for j=1:totalGenres
            % total Songs estimated to belong to particular genre
            numSongsEstAsGivenGenre = sum(Est_SongsIndxOfGivenGenre == j);
            confusionMatrix(genreID,j) = numSongsEstAsGivenGenre/totalSongsOfGenre;
        end
    end
    
    allConfusionMatrices(:,:,iter) = confusionMatrix;
    
end


% Computing Mean and Std-Dev of the Confusion Matrices
mean_ConfusionMatrix = 1/totalIter * sum(allConfusionMatrices,3)
% sum(mean_ConfusionMatrix,2)

diffSum = zeros(size(mean_ConfusionMatrix));
for i = 1: totalIter
    diffSum = diffSum + (allConfusionMatrices(:,:,i) - mean_ConfusionMatrix).^2;
end

stdDev_ConfusionMatrix = (1/totalIter*diffSum).^(1/2)
% sum(stdDev_ConfusionMatrix,2)
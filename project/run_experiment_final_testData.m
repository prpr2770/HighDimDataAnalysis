% run experiment 2: 
% training and test data


% read the training data

clear all; close all; clc;

% -------------------------------------------------------------------------------
% Parameters

trainDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
testDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\testTracks\';

totalIter = 20;              % totalIterations for Training-Test crossValids.
generate_ITML_DistanceMatrix = 0;   % Set 1: To RUN!

dataType = 'MFCC';
EXPERIMENT_TYPE = 'Graph'; % 'tSNE' 'Graph' 'ITML_Graph' 'ITML_tSNE'  'Original_Data'


% -------------------------------------------------------------------------------
% Read the training data

% matfile: storing 'CODEWORD_HIST_ALLSONGS'
aggregateDataDirName = fullfile(trainDirName,'aggregate_song_data')

% obtain cw-hist for all songs
cwHist_fileName = strcat(aggregateDataDirName,'\',dataType,'\cw_histogram_allSongs.mat');
cwHist_mFile = matfile(cwHist_fileName);
TRAIN_DATA = cwHist_mFile.CODEWORD_HIST_ALLSONGS';     % ROW-VECT DATASETS

% obtain songGenreID for trainingSongs
[genreKeys TRAIN_SONG_GENRES] = getGenreKeysForSongs(trainDirName);     % ROW-VEC

% obtain the DIST_METRIC
distMatrix_fileName = strcat(aggregateDataDirName,'\',dataType,'\itml_distance_matrix.mat');
distMatrix_mat = matfile(distMatrix_fileName, 'Writable', true);
DIST_METRIC = distMatrix_mat.DIST_METRIC;


% -------------------------------------------------------------------------------
% Read the testing data

aggregateDataDirName = fullfile(testDirName,'aggregate_song_data')
cwHist_fileName = strcat(aggregateDataDirName,'\cw_histogram_allSongs.mat');
cwHist_mFile = matfile(cwHist_fileName, 'Writable',true);
TEST_DATA = cwHist_mFile.CODEWORD_HIST_ALLSONGS';     % ROW-VECT DATASETS
% cwHist_mFile.CODEWORD_HIST_TRAINSONGS =  TRAIN_DATA;    % store the training_data into same file.

% obtain songGenreID for trainingSongs
songGenre_fileName = strcat(aggregateDataDirName,'\songGenres.mat');
songGenre_mFile = matfile(songGenre_fileName);
TEST_SONG_GENRES = songGenre_mFile.TESTSONG_GENRES;    %ROW-VEC

% create file to store the itml_distance_matrix
distMatrix_fileName = strcat(aggregateDataDirName,'\',dataType,'\itml_distance_matrix.mat');
distMatrix_mat = matfile(distMatrix_fileName, 'Writable', true);


% -------------------------------------------------------------------------------

ALL_DATA = [TRAIN_DATA; TEST_DATA];
ALL_SONG_GENRES = [TRAIN_SONG_GENRES TEST_SONG_GENRES];

size(ALL_DATA)
size(ALL_SONG_GENRES)

%{
%%  ITML :USE THE DIST_METRIC of TRAINING DATA
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Implement Experiments
% -------------------------------------------------------------------------
switch EXPERIMENT_TYPE
    
       
    case 'tSNE'
        % implement tSNE
        perplexity = 30;
        out_dims = 3;
        initial_dims = 30;
        
        DATA = tsne(ALL_DATA, ALL_SONG_GENRES, out_dims, initial_dims, perplexity);
        PLOT_TITLE = sprintf('tSNE Embedding of CodeWords-%s',dataType);
        
    case 'Graph'
        dimReducedSpace = 3;
        DATA = getRefinedGraphEmbedding(ALL_DATA, dimReducedSpace);
        PLOT_TITLE = sprintf('Refined Graph Embedding of CodeWords-%s',dataType);
        
    case 'ITML_tSNE'
        % read ITML_distance_matrix
        DIST_MATRIX = getDistanceMatrixFromMetric(ALL_DATA,DIST_METRIC);
        
        % implement tSNE
        perplexity = 30;
        out_dims = 3;
        
        DATA = tsne(DIST_MATRIX, ALL_SONG_GENRES, out_dims,  perplexity);
        PLOT_TITLE = sprintf('ITML + tSNE Embedding of CodeWords-%s',dataType);
        
    case 'ITML_Graph'
        % read ITML_distance_matrix
        DIST_MATRIX = getDistanceMatrixFromMetric(ALL_DATA,DIST_METRIC);       
       
        dimReducedSpace = 3;
        DATA = refinedGraphEmbedding_distanceMatrix(DIST_MATRIX, dimReducedSpace);
        
        PLOT_TITLE = sprintf('ITML + Graph Embedding of CodeWords-%s',dataType);
end

% scatter plot for all the data
fig11 = figure(11)
colormap jet
for i = 1:length(genreKeys)
    songsOfAGenre = find(songGenres == i);
    hue = (length(genreKeys)+1-i)*ones(length(songsOfAGenre),1);
    rad = 30*ones(length(songsOfAGenre),1);
    h = scatter3(DATA(songsOfAGenre,1),DATA(songsOfAGenre,2),DATA(songsOfAGenre,3),rad,hue,'filled');
    hold on
end
legend(genreKeys)
title(PLOT_TITLE)
hold off

figName = strcat(aggregateDataDirName,'\',dataType,'_',EXPERIMENT_TYPE,'_data.fig');
savefig(figName);
% %}

% ------------------------------------------------------------------------
%% Implement Clustering
% -------------------------------------------------------------------------


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

% detetmine length of training data
numTrainingSongs = length(TRAIN_DATA);

% ensure if it is column or row - matrices!     +++????????????????????????
anTrainingData = DATA(1:numTrainingSongs,:);
anTestData = DATA(numTrainingSongs+1:end,:);

anTestGenreIndex = reshape(TEST_SONG_GENRES,1,[]);
anTrainingGenreIndex = reshape(TRAIN_SONG_GENRES,1,[]);

BAM_test_experiment_results



% confusionMatrix = anTestData_ConfusionMatrix(:,:,4);
sumConf = sum(anConfusionMatrix,1);
factorMatrix = repmat(sumConf,length(anConfusionMatrix),1);
confusionMatrix = anConfusionMatrix./factorMatrix;

anPercentCorrect            % gives percent correct for the trainingSongs

nAllPercentCorrect          % gives allPercent correct for trainingSongs.

% ========================================================================




%}
%{
Experiments: CodeWord Histograms
%}
clear all; close all; clc;

% -------------------------------------------------------------------------------
% Parameters

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';

totalIter = 20;              % totalIterations for Training-Test crossValids.
generate_ITML_DistanceMatrix = 0;   % Set 1: To RUN!

dataType = 'MFCC';
EXPERIMENT_TYPE = 'Graph'; % 'tSNE' 'Graph' 'ITML_Graph' 'ITML_tSNE'  'Original_Data'

% -------------------------------------------------------------------------------
% FOLDER

% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'song_data')

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'aggregate_song_data')

% -------------------------------------------------------------------------------
% FILES

% matfile to store aggData_allSongs
aggData_allSongs_fileName = strcat(aggregateDataDirName,'\aggData_allSongs.mat');
aggData_mFile = matfile(aggData_allSongs_fileName,'Writable',true);

% matfile to store aggData_allSongs
results_fileName = strcat(aggregateDataDirName,'\results.mat');
results_mFile = matfile(results_fileName,'Writable',true);
% -------------------------------------------------------------------------

% matfile: storing 'CODEWORD_HIST_ALLSONGS'
cwHist_fileName = strcat(aggregateDataDirName,'\',dataType,'\cw_histogram_allSongs.mat');
cwHist_mFile = matfile(cwHist_fileName,'Writable',true);

distMatrix_fileName = strcat(aggregateDataDirName,'\',dataType,'\itml_distance_matrix.mat');
distMatrix_mat = matfile(distMatrix_fileName, 'Writable', true);

% -------------------------------------------------------------------------
% Obtain HighDimensional Dataset - Codeword Histogram

dataSet = cwHist_mFile.CODEWORD_HIST_ALLSONGS';     % ROW-VECT DATASETS

% obtain GenreIDs from TRUTH file.
[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

% -------------------------------------------------------------------------
if generate_ITML_DistanceMatrix
    % ---------------------------------------------------------------------
    % Obtain Distance_Metric and Distance_Matrix
    % ---------------------------------------------------------------------
    % implement iTML
    
    X = dataSet;       % ROW-MATRIX : Every row is a dataPoint
    y = songGenres';    % COL-VECTOR : Every row is a genreID for song.
    
    %     [dist_metric dist_matrix] = runITMLonDataSet(X,y)
    dist_metric = MetricLearningAutotuneKnn(@ItmlAlg, y, X);
    imagesc(dist_metric)
    colorbar
    distMatrix_mat.DIST_METRIC = dist_metric;
    
    % ------------------------------------------------------
    % compute Symmetric Distance_Matrix
    dist_matrix = getDistanceMatrixFromMetric(X,dist_metric);
    % ------------------------------------------------------
    % archive/save dist_matrix/metric
    distMatrix_mat.DIST_MATRIX = dist_matrix;
end

% %{
% =========================================================================
%% Experiments:
% -------------------------------------------------------------------------

switch EXPERIMENT_TYPE
    
    case 'Original_Data'
        DATA = dataSet;
        PLOT_TITLE = sprintf('Plot of CodeWords-%s',dataType);
        
    case 'tSNE'
        % implement tSNE
        perplexity = 30;
        out_dims = 3;
        initial_dims = 30;
        
        DATA = tsne(dataSet, songGenres, out_dims, initial_dims, perplexity);
        PLOT_TITLE = sprintf('tSNE Embedding of CodeWords-%s',dataType);
        
    case 'Graph'
        dimReducedSpace = 3;
        DATA = getRefinedGraphEmbedding(dataSet, dimReducedSpace);
        PLOT_TITLE = sprintf('Refined Graph Embedding of CodeWords-%s',dataType);
        
    case 'ITML_tSNE'
        % read ITML_distance_matrix
        dist_matrix = distMatrix_mat.DIST_MATRIX;
        
        % implement tSNE
        perplexity = 30;
        out_dims = 3;
        
        DATA = tsne(dist_matrix, songGenres, out_dims,  perplexity);
        PLOT_TITLE = sprintf('ITML + tSNE Embedding of CodeWords-%s',dataType);
        
    case 'ITML_Graph'
        % read ITML_distance_matrix
        dist_matrix = distMatrix_mat.DIST_MATRIX;
        
        % generate EMBEDDED points using GraphEmbedding
        %         DATA = getRefinedGraphEmbedding(Dataset, dimReducedSpace);
        
        dimReducedSpace = 3;
        DATA = refinedGraphEmbedding_distanceMatrix(dist_matrix, dimReducedSpace);
        
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


%% Implement Clustering
% -------------------------------------------------------------------------

[MEAN_C_MATRIX, STD_DEV_C_MATRIX, PROB_SUCCESS_MEAN, AVG_SUCCESS] =  run_experiment_train_test_BAM(tracksDirName, DATA, totalIter);
 

fig = figure()
colormap parula
subplot(2,1,1)
imagesc(MEAN_C_MATRIX);
title('MEAN Confusion Matrix')
colorbar

subplot(2,1,2)
imagesc(STD_DEV_C_MATRIX);
title('STD-DEV Confusion Matrix')
colorbar

figName = strcat(aggregateDataDirName,'\',dataType,'_',EXPERIMENT_TYPE,'_ConfMatrix.jpg');
saveas(fig,figName);

% -------------------------------------------------------------------------

msg = sprintf('Average Success: %s_CodeWords = %f',dataType,AVG_SUCCESS);
disp(msg);

eval(['results_mFile.',dataType,'_',EXPERIMENT_TYPE,'_PROB_SUCCESS_MEAN = PROB_SUCCESS_MEAN;'])
eval(['results_mFile.',dataType,'_',EXPERIMENT_TYPE,'_AVG_SUCCESS = AVG_SUCCESS;'])
eval(['results_mFile.',dataType,'_',EXPERIMENT_TYPE,'_STD_DEV_C_MATRIX = STD_DEV_C_MATRIX;'])
eval(['results_mFile.',dataType,'_',EXPERIMENT_TYPE,'_MEAN_C_MATRIX = MEAN_C_MATRIX;'])

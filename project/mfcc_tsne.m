% Initial Data Visualization and Processing
% 
% 1. Read and load the X datafile. 
% 2. appply visualziation of data
% 3. determine intrinsic dimension

clear all;close all; clc

dirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
matDataDir = fullfile(dirName,'matData_cut');

% Read the data files
mfcc_all_fileName = strcat(matDataDir,'\allSongs_mat_mfcc.mat');
mfcc_data = matfile(mfcc_all_fileName);

% load the X var
X = mfcc_data.X;        % true-Data colVectors
Y = mfcc_data.Y;        % true-labels (row-vec)



%{
% ======================================================================
% ----------------------------------------------------------------------
%  1. TSNE - Visualization Implementation

%   mappedX = tsne(X, labels, no_dims, initial_dims, perplexity)
%   mappedX = tsne(X, labels, initial_solution, perplexity)
%
% The function performs symmetric t-SNE on the NxD dataset X to reduce its 
% dimensionality to no_dims dimensions (default = 2).
% generate 


mappedX = tsne(X,[],3);

% visualization
dim1 = mappedX(:,1);
dim2 = mappedX(:,2);
dim3 = mappedX(:,3);

fig1 = figure(1)
plot3(dim1, dim2, dim3)

%}


% ======================================================================
%  2. INTRINSIC DIMENSION DETECTION

%%{
% Read the data files
mfcc_genre_mat_fileName = strcat(matDataDir,'\allGenre_Songs_mat_mfcc.mat');
mfccGenreFile = matfile(mfcc_genre_mat_fileName);

genreKeys = {'rock_pop', 'classical', 'electronic', 'jazz_blues', 'metal_punk', 'world'};

%%{

dimsGenre = [];

for genreIndex = 1:length(genreKeys)
    % read the genreDataFiles
    eval(['Genre' num2str(genreIndex) '= mfccGenreFile.Genre' num2str(genreIndex) ';'])
end


for genreIndex = 1:length(genreKeys)
    % Implement correlation dimension: :
    % GetDim(X) one column corresponds to one data point
    eval(['dim = GetDim(Genre' num2str(genreIndex) ');']);
    size(dim)
    dimsGenre = [dimsGenre; dim];
end

genreDimsDict = containers.Map(genreKeys, dimsGenre(:,2))


fig2 = figure(2)
plot(dimsGenre(:,1),'r')
hold on
plot(dimsGenre(:,2),'b')
hold on
plot(dimsGenre(:,3),'g')
hold off
title('Intrinsic Dimension of Datsets')

%%}

%{ 
Solution:

'classical'     [13.6232]
'electronic'    [7.5399]
'jazz_blues'    [16.0219]
'metal_punk'    [5.2997]
'rock_pop'      [6.4873]
'world'         [7.2136]

%}

%{
% ======================================================================
%  3. LEARN METRIC
% Create constraint set of the data: Similarity and Dissimilarity

disp('Starting ITML');
X = X';         % row-vector dataset
y = Y';         % truth labels: numeric labels

% size(X)
% size(y)

num_folds = 2;
knn_neighbor_size = 5;
tic
acc = CrossValidateKNN(y, X, @(y,X) MetricLearningAutotuneKnn(@ItmlAlg, y, X), num_folds, knn_neighbor_size);
toc
disp(sprintf('kNN cross-validated accuracy = %f', acc));

%}
% ======================================================================
% 1. mfcc_frames for different songs
% 2. GMM30 representation for each song
% 3. determine distance between the songs : loglikelihood
% 4. Build a similarity matrix between songs: Gaussian RBF
% 5. compute laplacian
% 6. determine eigenvectors

% extract dataset
dirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
matDataDir = fullfile(dirName,'matData_cut');
mfcc_all_fileName = strcat(matDataDir,'\class_artist_1_album_1_track_1_full_mfcc.mat');
mfcc_data = matfile(mfcc_all_fileName);

x = mfcc_data.MFCC;
x = x'; % row-vector dataset

% ==================================================
clear all; close all; clc
% Generate Data to Test
dims = 10; numCenters = 4;
% Generate the data
randn('state', 0); rand('state', 0);
gmix = gmm(dims, numCenters, 'spherical');
ndat1 = 20; ndat2 = 20; ndata = ndat1+ndat2;
gmix.centres =  repmat(rand(numCenters,1),1,dims)
gmix.covars = 0.1*ones(1,numCenters);
x = gmmsamp(gmix, ndata);

% ------------------------------------------
x_train = x(1:end-5,:);
x_test = x(end-4:end, :);

% ==================================================
% ----------------------------------------------
% Algo1: GMM-EM for clustering data
[numFrames numDims] = size(x_train)

% Set up mixture model
ncentres = 4; input_dim = numDims;
gmix = gmm(input_dim, ncentres, 'spherical')

% Initialise the mixture model
options = foptions;
gmix = gmminit(gmix, x_train, options)

% Loop over EM iterations.
% Expectation 
post = gmmpost(gmix, x_train)
% Maximization
options = foptions;
options(14) = 100; % Number of Iterations
options(1) = -1; % Switch off all messages, including warning
[gmix, options, errlog] = gmmem(gmix, x_train, options)


%GMMACTIV Computes the activations of a Gaussian mixture model.
activ_gmm = gmmactiv(gmix, x_test)
% activ_gmm = gmmactiv(gmix, x_train)
size(activ_gmm)
imagesc(activ_gmm)
colorbar

%GMMPROB Computes the data probability for a Gaussian mixture model.
prob_gmm = gmmprob(gmix, x_test(1,:))

%GMMPAK	Combines all the parameters in a Gaussian mixture model into one vector.
parameters_gmm = gmmpak(gmix)

% --------------------------------------------------
% Algo2: kNN+SingleGaussian_per_cluster

% k-means clustering - to detect clusters. Then fit a Single Gaussian to each cluster!
options = foptions;

[numSamples numDims] = size(x_train);
totalClusters = 5;
maxValues = 2 + ceil(max(x_train,[],1)); %column-wise maximum value to determine range of the centers
maxV_mat = repmat(maxValues,totalClusters,1);
minValues = -2 + floor(min(x_train,[],1));
minV_mat = repmat(minValues, totalClusters, 1);
centres = minV_mat + (maxV_mat - minV_mat).*rand(totalClusters, numDims); 

% run kmeans algo.
[centres, options, post] = kmeans(centres, x_train, options);

% set priors depending on points in each cluser
% Set priors depending on number of points in each cluster
cluster_sizes = max(sum(post, 1), 1);  % Make sure that no prior is zero
cluster_priors = cluster_sizes/sum(cluster_sizes); % Normalise priors

% how to fit SingleGaussian to a dataset?
% Determine MEAN, VARIANCE for each cluster
mean = zeros(totalClusters,numDims);
covarMatrix = zeros(numDims,numDims,totalClusters);
clusterPar = zeros(totalClusters,numDims  + numDims*numDims);
clusterWeights = cluster_priors';
for i = 1: totalClusters
    cluster = find(post(:,i));
    numPoints = length(cluster);
    if numPoints >0
        data = x_train(cluster,:);
        meanData = sum(data,1)/numPoints % colsum / numRows
        meanMat = repmat(meanData ,numPoints,1);
        covData = cov(data);
        covarVec = reshape(covData,1,[]);
        mean(i,:) = meanData;
        covarMatrix(:,:,i) = covData;
        clusterPar(i,:) = [meanData covarVec];
    else
        
    end
end

% Store the GMM representation of each song as .mat file!
% implement elementary code that does this computation. 


f1 = clusterPar
w1 = clusterWeights

% [counts binCenters] : histogram of two vectors
% [prior clusterCenter] : feature of a song
% f1, f2 -> [[cluserCenterMean Coords] clusterVar; ... ]
% w1, w2 -> priors
[f, fval] = emd(f1, f2, w1, w2, @gdf);

% mircluster and mirdistance
% c = mircluster(x_train)


% --------------------------------------------------
% Simple features: 
% Distance are computed as ABSOLUTE DIFFERENCES. L1 NORM.

% Computing distances between songs: when modelled by GMM 
% Monte-Carlo approach: Extract samples(N = Distance-Sample Rate) from the songA and songB
% determine the probability of these frame-samples being generated from GMM_SongB.
% Define a symmetric distance-measure:


% Distance between songs : EarthMoversDistance of the parameters of the GMM,
% KL-divergence between the GMMs representing the songs!
% Refer: Logan, Saloman, 2001 : A music similarity function based on signal analysis

% Comparison of Distance Measures between GMM
% https://www.ee.columbia.edu/~dpwe/pubs/JenECJ07-gmmdist.pdf



% Foote : Content based retrieval of Music and Audio
% mfcc -> Vector Quantization, Tree Quantization

% k-Means Feature Aggregation: Logan Saloman :http://www.hpl.hp.com/techreports/Compaq-DEC/CRL-2001-2.pdf

% Content based music similarity function :Beth Logan, Saloman

% ==================================================================
% For each frame, obtain a spectral representation.: mfcc
% 
% Cluster represetnation of song: 
% Each cluster: Mean, Variance, Weight
% {p_i} songA; {q_j} songB
% 
% Each cluster: {Mean, Sigma(covariance matrix), Weight of cluster}
% 
% d_pi_qj = distance between clusters pi and qj
% symmetric form of KL divergence: 
% d_pi_qj = S_pi / S_qj + S_qj/ S_pi + (mu_pi - mu_qj)^2 (1/S_pi + 1/S_qj)
% 
% 
% THus, we obtain a kxk matrix representing the distance between the songs. 
% f_pi_qj : kxk matrix -> Identify FLOW MATRIX that minimizes the total sum. 
% formulate as a LInear Programming. 
% 
% Compute EMD = ||D*F|| / ||F||
% ==================================================================

Fisher-VEctor Representation:
Refere: http://cmmr2012.eecs.qmul.ac.uk/sites/cmmr2012.eecs.qmul.ac.uk/files/pdf/papers/cmmr2012_submission_54.pdf




% ======================================================================
%  4. DIMENSION REDUCTION



% ======================================================================
%  5. KNN CLUSTERING ALGO


% ======================================================================
%  6. CLUSTER MEMBERSHIP


% Initial Data Visualization and Processing
% 
% 1. Read and load the X datafile. 
% 2. appply visualziation of data
% 3. determine intrinsic dimension


dirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
matDataDir = fullfile(dirName,'matData_cut');
status = mkdir(matDataDir);
mfcc_all_fileName = strcat(matDataDir,'\allSongs_mfcc.mat');
mfcc_data = matfile(mfcc_all_fileName);

% load the X var
X = mfcc_data.X;

% ======================================================================
% ----------------------------------------------------------------------
%  1. TSNE - Visualization Implementation

%   mappedX = tsne(X, labels, no_dims, initial_dims, perplexity)
%   mappedX = tsne(X, labels, initial_solution, perplexity)
%
% The function performs symmetric t-SNE on the NxD dataset X to reduce its 
% dimensionality to no_dims dimensions (default = 2).
% generate 


%{
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

genreKeys = {'rock_pop', 'classical', 'electronic', 'jazz_blues', 'metal_punk', 'world'};

mfcc_genre_fileName = strcat(matDataDir,'\allGenre_Songs_mfcc.mat');
mfccGenreFile = matfile(mfcc_genre_fileName);

dimsGenre = zeros(length(genreKeys),1);

for genreIndex = 1:length(genreKeys)
    % read the genreDataFiles
    eval(['Genre' num2str(genreIndex) '= mfccGenreFile.Genre' num2str(genreIndex) ';'])
end


for genreIndex = 1:length(genreKeys)
    % Implement correlation dimension: :
    % GetDim(X) one column corresponds to one data point
    eval(['dimsGenre(genreIndex) = GetDim(Genre' num2str(genreIndex) ');']);
end




[counts, centers] = hist(dimsGenre)

fig2 = figure(2)
subplot(2,1,1)
plot(centers, counts)
subplot(2,1,2)
plot(dimsGenre)

% ======================================================================
%  3. LEARN METRIC
% Create constraint set of the data: Similarity and Dissimilarity



% ======================================================================
%  4.  DIMENSION REDUCTION

% ======================================================================
%  5.  KNN CLUSTERING ALGO
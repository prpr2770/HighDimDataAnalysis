%{
Script to do the following:
1. Read the all_songs_dyn_mfcc to represent bag-of-frames: These are the
normalized dyn_mfcc_matrix
5. Extract the code-words from this bag-of-frames. ( k-means: 128,256, 512,...)
6. Save the code-word identifiers

7. Obtain code-word-histogram for each song.  (top-tau vector
quantization).
8. Store the code-word-histogram Identifiers of each song. 
%}

% =========================================================================
clear all; close all; clc;

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);


% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'dynMfcc_data')

% data directory to store mfcc and dyn_mfcc of all song
aggregateDataDirName = fullfile(tracksDirName,'dynMfcc_data_aggregate')

% read file
mfcc_allSongs_fileName = strcat(aggregateDataDirName,'\dyn_mfcc_allSongs.mat');
data_allSongs = load(mfcc_allSongs_fileName);

% create file: store norm_codeWordHistograms
cw_hist_allSongs_fileName = strcat(aggregateDataDirName,'\cw_hist_allSongs.mat');
cw_mFile = matfile(cw_hist_allSongs_fileName,'Writable',true);

% extract the global paramter/data
norm_dyn_mfcc_allSongs = data_allSongs.norm_dyn_mfcc_allSongs;
mean_dyn_mfcc = data_allSongs.mean_dyn_mfcc;
std_dyn_mfcc = data_allSongs.std_dyn_mfcc;
totalSongs = data_allSongs.countSong;

% =========================================================================
% Extract the codewords for all Songs:

numCodeWords  = 256;    % total CodeWord Clusters
tau = 10;               % numNearestNbrs of frame among CodeWordClusterCenters

[numCoeffs numFrames] = size(norm_dyn_mfcc_allSongs);

% Initialize randomly chosen nodes for kmeans
iter = 5;
X = norm_dyn_mfcc_allSongs';    % ROWVEC - Convert for ma_kmeans!
size(X)
[C,Qe,N,W,Q] = ma_kmeans(X, iter, numCodeWords);
codeWordCenters = C';           % COLVEC - Convert!

cw_mFile.codeWordCenters = codeWordCenters; % save(cw_hist_allSongs_fileName,'codeWordCenters');

% =========================================================================
% Extract the codewords for all Songs:

% read mfcc of each song from directory
% read all songs in the tracks folder
Files=dir(tracksDirName);
length(Files)



% norm_hist_codeWords_allSongs = zeros(numCodeWords,totalSongs);



countSong =0;
for k=1:length(Files)
    if (Files(k).bytes >0)
        countSong = countSong + 1
        wavfileName=Files(k).name;
        
        
        % Create file to save the mfcc and dyn_mfcc
        wavName = wavfileName(1:end-4);     % remove .wav extension
        mfcc_song_fileName = strcat(matDataDirName,'\',wavName,'.mat');
        
        mfcc_song_File = matfile(mfcc_song_fileName);
        NORM_DYN_MFCC = mfcc_song_File.NORM_DYN_MFCC;
        
        %{
        % --------------------------------------------
        % normalize the dyn_mfcc_song
        [numCoeffs numFrames] = size(dyn_mfcc_song);
        
        mean_dyn_mfcc_mat = repmat(mean_dyn_mfcc,1,numFrames);
        std_dyn_mfcc_mat = repmat(std_dyn_mfcc,1,numFrames);
        
        norm_dyn_mfcc_song = (dyn_mfcc_song - mean_dyn_mfcc_mat)./std_dyn_mfcc_mat;
        % --------------------------------------------
        %}


        % --------------------------------------------
        % determine the tau-nearest-nbrs for each frame-in-song among
        % codeWordCenters.
        size(codeWordCenters)
        size(NORM_DYN_MFCC)
        nbrs_of_songFrames = knnsearch(codeWordCenters',NORM_DYN_MFCC','k',tau,'distance','euclidean');
        nbrs_of_songFrames = reshape(nbrs_of_songFrames,[],1);
        
        % extract histogram of occurence and re-structure into COL_VEC
        hist_vec = histc(nbrs_of_songFrames,1:numCodeWords);
        hist_vec = reshape(hist_vec,[],1);
        
        if countSong > 1
        % update the norm_hist_codeWords representation of each song.
        % norm_hist_codeWords_allSongs(:,countSong) = (1/numFrames) * (1/tau)* hist_vec;

        % save directly into file
        cw_mFile.NORM_CW_HISTS(:,countSong) = (1/numFrames) * (1/tau)* hist_vec;
    else
        % save into matrix
        % norm_hist_codeWords_allSongs(:,countSong) = (1/numFrames) * (1/tau)* hist_vec;

        % save directly into file
        cw_mFile.NORM_CW_HISTS = (1/numFrames) * (1/tau)* hist_vec;
        end


        % iteratively save the norm_hist_codeWord into the file.
        
    end % iterate over all songs
end %iterate over all files

% archive the normalized_codeWord_histogram_representation.
% save(cw_hist_allSongs_fileName,'norm_hist_codeWords_allSongs','-append');

size(cw_mFile.NORM_CW_HISTS)
imagesc(cw_mFile.NORM_CW_HISTS)
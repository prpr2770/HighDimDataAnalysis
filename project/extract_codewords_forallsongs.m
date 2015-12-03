%{
Script that does the following:
+ Read CODEWORDS_BATCH_ALLSONGS
+ Extract the final CODEWORDS from this batch
+ Save the CODEWORDS in a separate file.


%}


function extract_codewords_forallsongs(numCodeWords, iter)

% numCodeWords = 256;
% iter = 100;
% tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';

   
   [genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);
    
    % data directory to store mfcc and dyn_mfcc of each song
    matDataDirName = fullfile(tracksDirName,'song_data')
    
    % data directory to store global/aggregate information
    aggregateDataDirName = fullfile(tracksDirName,'aggregate_song_data')
    
    
    % -------------------------------------------------------------------------------
    
    % matfile to store aggData_allSongs
    aggData_allSongs_fileName = strcat(aggregateDataDirName,'\aggData_allSongs.mat');
    aggData_mFile = matfile(aggData_allSongs_fileName,'Writable',true);
    
    % matfile: storing 'CODEWORDS_BATCH_ALLSONGS'
    codewords_fileName = strcat(aggregateDataDirName,'\codewords_batch_allSongs.mat');
    codewords_mFile = matfile(codewords_fileName,'Writable',true);


% matfile: storing 'CODEWORDS_ALLSONGS'
final_codewords_fileName = strcat(aggregateDataDirName,'\codewords_allSongs.mat');
final_codewords_mFile = matfile(final_codewords_fileName,'Writable',true);

% =============================================================================================

% read all the codewords_batch_allsongs.
CODEWORDS_BATCH_ALLSONGS = codewords_mFile.CODEWORDS_BATCH_ALLSONGS;

% -------------------------------------------------------------------------------
% Implement the kmeans algorithm

CODEWORDS_BATCH_ALLSONGS = CODEWORDS_BATCH_ALLSONGS';		% row-vec Dataset

tic
CODEWORDS_ALLSONGS = ma_kmeans(CODEWORDS_BATCH_ALLSONGS, iter, numCodeWords);
toc

CODEWORDS_ALLSONGS = CODEWORDS_ALLSONGS';					% COL-VEC Dataset

% -------------------------------------------------------------------------------
% store the final_codewords
final_codewords_mFile.CODEWORDS_ALLSONGS = CODEWORDS_ALLSONGS;

clear all; close all; clc;
end
% =============================================================================================
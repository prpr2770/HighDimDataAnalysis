%{
Script that does the following:
+ Reads all the CODEWORDS
+ Reads every individual song, and obtain the cw-historgram representation of each song.
+ Save the code-words of each song in a different-file.

%}


function extract_codeWord_Histograms( numCodeWords, tau)
% =========================================================================
% numCodeWords = 256;
% tau = 10;

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'dynMfcc_data')

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'dynMfcc_data_aggregate')


% -------------------------------------------------------------------------------
% matfile to store aggData_allSongs
aggData_allSongs_fileName = strcat(aggregateDataDirName,'\aggData_allSongs.mat');
aggData_mFile = matfile(aggData_allSongs_fileName);

% matfile: storing 'CODEWORDS_ALLSONGS'
codewords_fileName = strcat(aggregateDataDirName,'\codewords_allSongs.mat');
codewords_mFile = matfile(codewords_fileName);

% matfile: storing 'CODEWORD_HIST_ALLSONGS'
cwHist_fileName = strcat(aggregateDataDirName,'\cw_histogram_allSongs.mat');
cwHist_mFile = matfile(cwHist_fileName,'Writable',true);


% =============================================================================================
% Read all the codewords_mFile
CODEWORDS_ALLSONGS = codewords_mFile.CODEWORDS_ALLSONGS;


% --------------------------------------------------------------
% Determine tau-nearest-neighbor CodeWord_Histogram for each song.

Files = dir(matDataDirName);

countSong = 0;
for k=1:length(Files)       % sequentially analyze dyn_mfcc_data song-wise
    
    if (Files(k).bytes > 0 )
        countSong = countSong + 1
        
        % Read the mat file containing mfcc,dyn_mfcc data
        matfileName=Files(k).name;
        mfcc_song_fileName = strcat(matDataDirName,'\',matfileName);
        mfcc_song_mFile = matfile(mfcc_song_fileName);
        
        
        % extract the data needed.
        NORM_DYN_MFCC = mfcc_song_mFile.NORM_DYN_MFCC;
        [coeffs songFrames] = size(NORM_DYN_MFCC);
        
        
        % find tau-nearest neighbors
        
        size(CODEWORDS_ALLSONGS)
        size(NORM_DYN_MFCC)
        nbrs_of_songFrames = knnsearch(CODEWORDS_ALLSONGS',NORM_DYN_MFCC','k',tau,'distance','euclidean');
        nbrs_of_songFrames = reshape(nbrs_of_songFrames,[],1);
        
        % extract histogram of occurence and re-structure into COL_VEC
        hist_vec = histc(nbrs_of_songFrames,1:numCodeWords);
        hist_vec = reshape(hist_vec,[],1);
        
        % save directly into file
        if countSong > 1
            cwHist_mFile.CODEWORD_HIST_ALLSONGS(:,countSong) = (1/songFrames) * (1/tau)* hist_vec;
        else
            cwHist_mFile.CODEWORD_HIST_ALLSONGS = (1/songFrames) * (1/tau)* hist_vec;
        end
        
        
        
    end
end

clear all; close all; clc;
end

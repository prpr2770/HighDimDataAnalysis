function extract_MFCC_CodeWords_BatchOfSongs(numCodeWordsPerBatch, songsPerBatch, iter)
%{
	Script that does the following:
	+ Randomly select 10 songs files's NORM_MFCC
	+ Implement k-means, to detect codeWords
	+ Archive CodeWord vectors into a GlobalRepository, for each collection of 10 songs


% numCodeWordsPerBatch = 1024;
% songsPerBatch = 10;
% iter = 50;              % kMeans algorithm

%}

% =========================================================================

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'dynMfcc_data')

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'dynMfcc_data_aggregate')


% -------------------------------------------------------------------------------

% matfile to store aggData_allSongs
aggData_allSongs_fileName = strcat(aggregateDataDirName,'\aggData_allSongs.mat');
aggData_mFile = matfile(aggData_allSongs_fileName,'Writable',true);

% 
% % matfile: storing 'NRM_DYN_MFCC_ALLSONGS'
% nrm_dyn_mfcc_allSongs_fileName = strcat(aggregateDataDirName,'\nrm_dyn_mfcc_allSongs.mat');
% nrm_dyn_mfcc_allSongs_mFile = matfile(nrm_dyn_mfcc_allSongs_fileName,'Writable',true);

% matfile: storing 'CODEWORDS_BATCH_ALLSONGS'
codewords_fileName = strcat(aggregateDataDirName,'\codewords_mfcc_batch_allSongs.mat');
codewords_mFile = matfile(codewords_fileName,'Writable',true);

% =============================================================================================
% =============================================================================================

% Obtain Aggregate Information
TOTAL_FRAMES = aggData_mFile.TOTAL_FRAMES;
TOTAL_SONGS = aggData_mFile.TOTAL_SONGS ;
DIMS_MFCC = aggData_mFile.dims_MFCC ;

% ---------------------------------------------------------------------------------------------
% Compute CodeWords for batches of songs

% numCodeWordsPerBatch = 1024;
% songsPerBatch = 10;
% iter = 50;              % kMeans algorithm

Files = dir(matDataDirName);
rndOrder_songs = randperm(length(Files));

totalBatches = ceil(length(Files)/songsPerBatch);

fileCount = 1;
batchCount = 1;
songCount = 0;

while (fileCount < length(Files) + 1)
    
    % ------------------------------------------------------------------------------
    % concatanate all norm_dyn_mfcc frames for each batch.
    
    BATCH_NORM_MFCC = [];
    
    % concatanate
    while(mod(fileCount, songsPerBatch) ~=0 && fileCount <= length(Files) )
        warning('Concatanating songs data')
        songID = rndOrder_songs(fileCount);
        
        if (Files(songID).bytes > 0 )
            songCount = songCount + 1;
            % Read the mat file containing mfcc,dyn_mfcc data
            matfileName=Files(songID).name;
            mfcc_song_fileName = strcat(matDataDirName,'\',matfileName);
            mfcc_song_mFile = matfile(mfcc_song_fileName);
            
            
            % extract the data needed.
            NORM_MFCC = mfcc_song_mFile.NORM_MFCC;
            
            % concatanate the frames for the 10 songs
            BATCH_NORM_MFCC = [BATCH_NORM_MFCC NORM_MFCC];
            
        end
        % increment fileCount
        fileCount = fileCount+1
    end
    
    % ------------------------------------------------------------------------------
    % Implement the kmeans algorithm to detect cluster
    size(BATCH_NORM_MFCC)
    tic
    BATCH_CODEWORDS = ma_kmeans(BATCH_NORM_MFCC', iter, numCodeWordsPerBatch);
    toc
    BATCH_CODEWORDS = BATCH_CODEWORDS';           % COLVEC - Convert!
    
    % ------------------------------------------------------------------------------
    % store/archive the clusters computed.
    % -------------------------------------------------------
    % archive the norm_dyn_mfcc_ofAllSongs!
    if batchCount > 1
        [nrows ncols] = size(codewords_mFile,'CODEWORDS_BATCH_ALLSONGS');
        numFrames = size(BATCH_CODEWORDS,2);
        codewords_mFile.CODEWORDS_BATCH_ALLSONGS(:,ncols+1:ncols+numFrames) = BATCH_CODEWORDS;
    else
        codewords_mFile.CODEWORDS_BATCH_ALLSONGS = BATCH_CODEWORDS;
    end
    
    
    
    % increment fileCount
    fileCount = fileCount+1
    batchCount = batchCount + 1
end


clear all; close all; clc;
end
%{
Script to do the following:
1. Read all the .mat files containing MFCC and DYN_MFCC and do the following:
a) Compute the MEAN_MFCC
b) Compute STD_MFCC
c) Compute and store NORM_MFCC in the same .mat file.

2. Create a globalData_store containing following information:
a) MEAN_MFCC
b) STD_MFCC



%}

% =========================================================================
clear all; close all; clc;

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



% =============================================================================================
% =============================================================================================

% Obtain Aggregate Information
TOTAL_FRAMES = aggData_mFile.TOTAL_FRAMES;
TOTAL_SONGS = aggData_mFile.TOTAL_SONGS ;
DIMS_DYN_MFCC = aggData_mFile.dims_DYN_MFCC ;
DIMS_MFCC = aggData_mFile.dims_DYN_MFCC ;

% -------------------------------------------------------------------------------
% Compute the MEAN_DYN_MFCC : song-wise
sumMean = zeros(DIMS_MFCC,1);

countSong = 0;
totalFrames = 0;            % track total frames evaluated, when analysing song-wise.

Files = dir(matDataDirName);
for k=1:length(Files)       % sequentially analyze dyn_mfcc_data song-wise
    if (Files(k).bytes >0 )
        countSong = countSong + 1
        
        % Read the mat file containing mfcc,dyn_mfcc data
        matfileName=Files(k).name;
        mfcc_song_fileName = strcat(matDataDirName,'\',matfileName);
        mfcc_song_mFile = matfile(mfcc_song_fileName);
        

        % extract the data needed.
        MFCC = mfcc_song_mFile.MFCC;
        [coeffs songFrames] = size(MFCC);
        
        sumMean = sum(MFCC,2);  %rowsum
        totalFrames = totalFrames + songFrames;
        
        % clear memory
        clear MFCC;        
    end
end
       
mean_mfcc = (1/totalFrames) * sumMean;
size(mean_mfcc)

MEAN_MFCC = mean_mfcc;

aggData_mFile.MEAN_MFCC = MEAN_MFCC;

% =============================================================================================
% =============================================================================================
% Compute STD_DYN_MFCC : song-wise

varSum = zeros(DIMS_MFCC,1);

countSong = 0;
totalFrames = 0;            % track total frames evaluated, when analysing song-wise.

Files = dir(matDataDirName);
for k=1:length(Files)       % sequentially analyze dyn_mfcc_data song-wise
    if (Files(k).bytes >0 )
        countSong = countSong + 1
        
        % Read the mat file containing mfcc,dyn_mfcc data
        matfileName=Files(k).name;
        mfcc_song_fileName = strcat(matDataDirName,'\',matfileName);
        mfcc_song_mFile = matfile(mfcc_song_fileName);
        

        % extract the data needed.
        MFCC = mfcc_song_mFile.MFCC;
        [coeffs songFrames] = size(MFCC);

        mean_mfcc_mat = repmat(mean_mfcc,1,songFrames);

        varSum = sum((MFCC - mean_mfcc_mat).^2,2);  %rowsum
        totalFrames = totalFrames + songFrames;

        clear MFCC mean_mfcc_mat;
        
    end
end
       
std_mfcc = ((1/totalFrames) * varSum).^(0.5);
size(std_mfcc)

STD_MFCC = std_dyn_mfcc;


aggData_mFile.STD_MFCC = STD_MFCC;


% %{
% ===================================================================================================
% computing and storing the normalized_dyn_mfcc

warning('Computing and storing Normalized Dyn_MFCC for ALL_FRAMES');


tic

countSong = 0;
totalFrames = 0;            % track total frames evaluated, when analysing song-wise.

Files = dir(matDataDirName);

for k=1:length(Files)       % sequentially analyze dyn_mfcc_data song-wise

    if (Files(k).bytes > 0 )
        countSong = countSong + 1
        
        % Read the mat file containing mfcc,dyn_mfcc data
        matfileName=Files(k).name;
        mfcc_song_fileName = strcat(matDataDirName,'\',matfileName);
        mfcc_song_mFile = matfile(mfcc_song_fileName,'Writable',true);
        

        % extract the data needed.
        MFCC = mfcc_song_mFile.MFCC;
        [coeffs songFrames] = size(MFCC);

        % -------------------------------------------------------
        % compute the NORM-DYN-MFCC
        mean_mfcc_mat = repmat(mean_mfcc,1,songFrames);
        std_mfcc_mat = repmat(std_dyn_mfcc,1,songFrames);
        NORM_MFCC = (MFCC - mean_mfcc_mat)./std_mfcc_mat;

        % store the NORM_DYN_MFCC
        mfcc_song_mFile.NORM_MFCC = NORM_MFCC;

        % -------------------------------------------------------
        % archive the norm_dyn_mfcc_ofAllSongs!
        if countSong > 1
            [nrows ncols] = size(nrm_dyn_mfcc_allSongs_mFile,'NRM_MFCC_ALLSONGS');
            numFrames = size(NORM_MFCC,2);
            nrm_dyn_mfcc_allSongs_mFile.NRM_MFCC_ALLSONGS(:,ncols+1:ncols+numFrames) = NORM_MFCC;
        else
            nrm_dyn_mfcc_allSongs_mFile.NRM_MFCC_ALLSONGS = NORM_MFCC;
        end

        % -------------------------------------------------------
        % clear up memory
        clear mean_mfcc_mat std_mfcc_mat NORM_MFCC MFCC ;
        
    end
end
       
toc


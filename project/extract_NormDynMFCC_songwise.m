%{
Script to do the following:
1. Read all the .mat files containing MFCC and DYN_MFCC and do the following:
a) Compute the MEAN_DYN_MFCC
b) Compute STD_DYN_MFCC
c) Compute and store NORM_DYN_MFCC in the same .mat file.

2. Create a globalData_store containing following information:
a) MEAN_DYN_MFCC
b) STD_DYN_MFCC
c) totalSongs


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


% matfile: storing 'NRM_DYN_MFCC_ALLSONGS'
nrm_dyn_mfcc_allSongs_fileName = strcat(aggregateDataDirName,'\nrm_dyn_mfcc_allSongs.mat');
nrm_dyn_mfcc_allSongs_mFile = matfile(nrm_dyn_mfcc_allSongs_fileName,'Writable',true);


% =============================================================================================
% =============================================================================================

% Obtain Aggregate Information
TOTAL_FRAMES = aggData_mFile.TOTAL_FRAMES;
TOTAL_SONGS = aggData_mFile.TOTAL_SONGS ;
DIMS_DYN_MFCC = aggData_mFile.dims_DYN_MFCC ;

% -------------------------------------------------------------------------------
% Compute the MEAN_DYN_MFCC : song-wise
sumMean = zeros(DIMS_DYN_MFCC,1);

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
        DYN_MFCC = mfcc_song_mFile.DYN_MFCC;
        [coeffs songFrames] = size(DYN_MFCC);
        
        sumMean = sum(DYN_MFCC,2);  %rowsum
        totalFrames = totalFrames + songFrames;
        
        % clear memory
        clear DYN_MFCC;        
    end
end
       
mean_dyn_mfcc = (1/totalFrames) * sumMean;
size(mean_dyn_mfcc)

MEAN_DYN_MFCC = mean_dyn_mfcc;

aggData_mFile.MEAN_DYN_MFCC = MEAN_DYN_MFCC;

% =============================================================================================
% =============================================================================================
% Compute STD_DYN_MFCC : song-wise

varSum = zeros(DIMS_DYN_MFCC,1);

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
        DYN_MFCC = mfcc_song_mFile.DYN_MFCC;
        [coeffs songFrames] = size(DYN_MFCC);

        mean_dyn_mfcc_mat = repmat(mean_dyn_mfcc,1,songFrames);

        varSum = sum((DYN_MFCC - mean_dyn_mfcc_mat).^2,2);  %rowsum
        totalFrames = totalFrames + songFrames;

        clear DYN_MFCC mean_dyn_mfcc_mat;
        
    end
end
       
std_dyn_mfcc = ((1/totalFrames) * varSum).^(0.5);
size(std_dyn_mfcc)

STD_DYN_MFCC = std_dyn_mfcc;


aggData_mFile.STD_DYN_MFCC = STD_DYN_MFCC;


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
        DYN_MFCC = mfcc_song_mFile.DYN_MFCC;
        [coeffs songFrames] = size(DYN_MFCC);

        % -------------------------------------------------------
        % compute the NORM-DYN-MFCC
        mean_dyn_mfcc_mat = repmat(mean_dyn_mfcc,1,songFrames);
        std_dyn_mfcc_mat = repmat(std_dyn_mfcc,1,songFrames);
        NORM_DYN_MFCC = (DYN_MFCC - mean_dyn_mfcc_mat)./std_dyn_mfcc_mat;

        % store the NORM_DYN_MFCC
        mfcc_song_mFile.NORM_DYN_MFCC = NORM_DYN_MFCC;

        % -------------------------------------------------------
        % archive the norm_dyn_mfcc_ofAllSongs!
        if countSong > 1
            [nrows ncols] = size(nrm_dyn_mfcc_allSongs_mFile,'NRM_DYN_MFCC_ALLSONGS');
            numFrames = size(NORM_DYN_MFCC,2);
            nrm_dyn_mfcc_allSongs_mFile.NRM_DYN_MFCC_ALLSONGS(:,ncols+1:ncols+numFrames) = NORM_DYN_MFCC;
        else
            nrm_dyn_mfcc_allSongs_mFile.NRM_DYN_MFCC_ALLSONGS = NORM_DYN_MFCC;
        end

        % -------------------------------------------------------
        % clear up memory
        clear mean_dyn_mfcc_mat std_dyn_mfcc_mat NORM_DYN_MFCC DYN_MFCC ;
        
    end
end
       
toc


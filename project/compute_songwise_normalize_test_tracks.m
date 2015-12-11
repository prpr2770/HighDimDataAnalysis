%{
Script to do the following:
1. Read all the .mat files containing MFCC, FD_MFCC and DYN_MFCC and do the following:
a) Compute the MEAN_DYN_MFCC
b) Compute STD_DYN_MFCC
c) Compute and store NORM_DYN_MFCC in the same .mat file.

2. Create a globalData_store containing following information:
a) MEAN_[]
b) STD_[]

%}

% =========================================================================
% clear all; close all; clc;
% tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';

% function compute_songwise_normalized_data(tracksDirName)


% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'song_data')

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'aggregate_song_data')


% -------------------------------------------------------------------------------
% matfile to store aggData_allSongs
aggData_allSongs_fileName = strcat(aggregateDataDirName,'\aggData_allSongs.mat');
aggData_mFile = matfile(aggData_allSongs_fileName,'Writable',true);

% -------------------------------------------------------------------------------
% Obtain Aggregate Information
TOTAL_FRAMES = aggData_mFile.TOTAL_FRAMES;
TOTAL_SONGS = aggData_mFile.TOTAL_SONGS ;

DIMS_DYN_MFCC = aggData_mFile.dims_DYN_MFCC ;
DIMS_MFCC = aggData_mFile.dims_MFCC ;
DIMS_FD_DYN_MFCC = aggData_mFile.dims_FD_DYN_MFCC ;

% =============================================================================================
% =============================================================================================


% Read the MEAN and STD_DEV from the Training dataset.
train_aggData_fileName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\aggregate_song_data\aggData_allSongs.mat';
train_aggData_mFile = matfile(train_aggData_fileName);

MEAN_MFCC = train_aggData_mFile.MEAN_MFCC ;
MEAN_DYN_MFCC = train_aggData_mFile.MEAN_DYN_MFCC ;
MEAN_FD_DYN_MFCC = train_aggData_mFile.MEAN_FD_DYN_MFCC;

STD_MFCC = train_aggData_mFile.STD_MFCC;
STD_DYN_MFCC = train_aggData_mFile.STD_DYN_MFCC ;
STD_FD_DYN_MFCC = train_aggData_mFile.STD_FD_DYN_MFCC;


% -------------------------------------------------------------------------------
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
        
       
        % ----------------------------------------------------
        % extract mfcc
        DATA = mfcc_song_mFile.MFCC;
        MEAN = MEAN_MFCC;
        STD = STD_MFCC;
        
        % compute the normalized data
        [coeffs songFrames] = size(DATA);
        MEAN_mat = repmat(MEAN,1,songFrames);
        STD_mat = repmat(STD,1,songFrames);
        NORM_DATA = (DATA - MEAN_mat)./STD_mat;

        % store the data
        mfcc_song_mFile.NORM_MFCC = NORM_DATA;

        
        % ----------------------------------------------------
        % extract the data needed.
        DATA = mfcc_song_mFile.DYN_MFCC;
        MEAN = MEAN_DYN_MFCC;
        STD = STD_DYN_MFCC;
        
        % compute the normalized data
        [coeffs songFrames] = size(DATA);
        MEAN_mat = repmat(MEAN,1,songFrames);
        STD_mat = repmat(STD,1,songFrames);
        NORM_DATA = (DATA - MEAN_mat)./STD_mat;

        % store the data
        mfcc_song_mFile.NORM_DYN_MFCC = NORM_DATA;

        % ----------------------------------------------------
        % extract the data needed.
        DATA = mfcc_song_mFile.FD_DYN_MFCC;
        MEAN = MEAN_FD_DYN_MFCC;
        STD = STD_FD_DYN_MFCC;
        
        % compute the normalized data
        [coeffs songFrames] = size(DATA);
        MEAN_mat = repmat(MEAN,1,songFrames);
        STD_mat = repmat(STD,1,songFrames);
        NORM_DATA = (DATA - MEAN_mat)./STD_mat;

        % store the data
        mfcc_song_mFile.NORM_FD_DYN_MFCC = NORM_DATA;

        
    end
end
       
toc

% end
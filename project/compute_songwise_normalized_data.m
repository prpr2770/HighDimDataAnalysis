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


% =============================================================================================
% =============================================================================================

% Obtain Aggregate Information
TOTAL_FRAMES = aggData_mFile.TOTAL_FRAMES;
TOTAL_SONGS = aggData_mFile.TOTAL_SONGS ;

DIMS_DYN_MFCC = aggData_mFile.dims_DYN_MFCC ;
DIMS_MFCC = aggData_mFile.dims_MFCC ;
DIMS_FD_DYN_MFCC = aggData_mFile.dims_FD_DYN_MFCC ;

% -------------------------------------------------------------------------------
% Compute the MEAN_DYN_MFCC : song-wise


sumMean_MFCC = zeros(DIMS_MFCC,1);
sumMean_FD_DYN_MFCC = zeros(DIMS_FD_DYN_MFCC,1);
sumMean_DYN_MFCC = zeros(DIMS_DYN_MFCC,1);

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
        
        % ----------------------------------------------------
        % extract the data needed.
        MFCC = mfcc_song_mFile.MFCC;
        [coeffs songFrames] = size(MFCC);
        sumMean_MFCC = sum(MFCC,2);  %rowsum
        
        % ----------------------------------------------------
        % extract the data needed.
        DYN_MFCC = mfcc_song_mFile.DYN_MFCC;
        [coeffs songFrames] = size(DYN_MFCC);
        sumMean_DYN_MFCC = sum(DYN_MFCC,2);  %rowsum
        
                
        % ----------------------------------------------------
        % extract the data needed.
        FD_DYN_MFCC = mfcc_song_mFile.FD_DYN_MFCC;
        [coeffs songFrames] = size(FD_DYN_MFCC);
        sumMean_FD_DYN_MFCC = sum(FD_DYN_MFCC,2);  %rowsum
        
        % ----------------------------------------------------
        totalFrames = totalFrames + songFrames;
        
        % clear memory
        clear DYN_MFCC MFCC FD_DYN_MFCC;        
    end
end
       

MEAN_MFCC = (1/totalFrames) * sumMean_MFCC;
MEAN_FD_DYN_MFCC = (1/totalFrames) * sumMean_FD_DYN_MFCC ;
MEAN_DYN_MFCC = (1/totalFrames) * sumMean_DYN_MFCC ;

aggData_mFile.MEAN_MFCC = MEAN_MFCC;
aggData_mFile.MEAN_DYN_MFCC = MEAN_DYN_MFCC;
aggData_mFile.MEAN_FD_DYN_MFCC = MEAN_FD_DYN_MFCC;
% =============================================================================================
% =============================================================================================
% Compute STD_DYN_MFCC : song-wise

varSum_MFCC = zeros(DIMS_MFCC,1);
varSum_FD_DYN_MFCC = zeros(DIMS_FD_DYN_MFCC,1);
varSum_DYN_MFCC = zeros(DIMS_DYN_MFCC,1);

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
        
        % ----------------------------------------------------
        % extract the data needed.
        MFCC = mfcc_song_mFile.MFCC;
        [coeffs songFrames] = size(MFCC);
        
        MEAN_MFCC_mat = repmat(MEAN_MFCC,1,songFrames);
        varSum_MFCC = sum((MFCC - MEAN_MFCC_mat).^2,2);  %rowsum

        % ----------------------------------------------------
        % extract the data needed.
        DYN_MFCC = mfcc_song_mFile.DYN_MFCC;
        [coeffs songFrames] = size(DYN_MFCC);
        
        MEAN_DYN_MFCC_mat = repmat(MEAN_DYN_MFCC,1,songFrames);
        varSum_DYN_MFCC = sum((DYN_MFCC - MEAN_DYN_MFCC_mat).^2,2);  %rowsum
                
        % ----------------------------------------------------
        % extract the data needed.
        FD_DYN_MFCC = mfcc_song_mFile.FD_DYN_MFCC;
        [coeffs songFrames] = size(FD_DYN_MFCC);

        MEAN_FD_DYN_MFCC_mat = repmat(MEAN_FD_DYN_MFCC,1,songFrames);
        varSum_FD_DYN_MFCC = sum((FD_DYN_MFCC - MEAN_FD_DYN_MFCC_mat).^2,2);  %rowsum

        % ----------------------------------------------------
        
        
        
        totalFrames = totalFrames + songFrames;

        clear DYN_MFCC MFCC FD_DYN_MFCC MEAN_FD_DYN_MFCC_mat MEAN_DYN_MFCC_mat MEAN_MFCC_mat;
        
    end
end


STD_MFCC = ((1/totalFrames) * varSum_MFCC).^(0.5);
STD_DYN_MFCC = ((1/totalFrames) * varSum_DYN_MFCC).^(0.5);
STD_FD_DYN_MFCC = ((1/totalFrames) * varSum_FD_DYN_MFCC).^(0.5);

aggData_mFile.STD_MFCC = STD_MFCC;
aggData_mFile.STD_DYN_MFCC = STD_DYN_MFCC;
aggData_mFile.STD_FD_DYN_MFCC = STD_FD_DYN_MFCC;


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
%{
Script to do the following:
Read all the songs in tracks.
Extract MFCC, FD_MFCC, DYN_MFCC information.
Store all these into file for a single-song.
%}

% =========================================================================
% clear all; close all; clc;


function extract_songwise_data(tracksDirName)

% tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';

[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);




% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'song_data')
status = mkdir(matDataDirName)

% data directory to store mfcc and dyn_mfcc of all song
aggregateDataDirName = fullfile(tracksDirName,'aggregate_song_data')
status = mkdir(aggregateDataDirName)


% ----------------------------------------------------------------------------------------------
% matfile to store aggData_allSongs
aggData_allSongs_fileName = strcat(aggregateDataDirName,'\aggData_allSongs.mat');
aggData_mFile = matfile(aggData_allSongs_fileName,'Writable',true);


% read all songs and extract mfcc
Files=dir(tracksDirName);
length(Files)

% ===============================================================================================================
% ===============================================================================================================
% Extract and store MFCC and DYN_MFCC data.
% -------------------------------------------------------------------------
% p = struct('fs',11025,'visu',0,'fft_size',256,'hopsize',128,'num_ceps_coeffs',0,'dB_max',96,'mel_filt_bank','auditory-toolbox');
p = struct('fs',11025,'visu',0,'fft_size',256,'hopsize',128,'num_ceps_coeffs',20,'dB_max',96);


countSong = 0;              % count total-songs processed in DB
totalFrames = 0;            % count total-frames generated for all songs
for k=1: length(Files)
    if (Files(k).bytes >0 )
        countSong = countSong + 1
        wavfileName=Files(k).name;
        
        
        % Create file to save the mfcc and dyn_mfcc
        wavName = wavfileName(1:end-4);     % remove .wav extension
        mfcc_song_fileName = strcat(matDataDirName,'\',wavName,'.mat');
        
        
        % -------------------------------------------------------
        % start reading the files
        
        % construct full file-name
        wavFullFileName = fullfile(tracksDirName,wavfileName);
        
        try
            % read audio-file
            [wavFile, Fs] = audioread(wavFullFileName);
            wavInfo = audioinfo(wavFullFileName);
            
        catch E1
            msg = strcat('Error in audioread: ',wavFullFileName);
            warning(msg)
        end
        
        % generate mfcc and dct
        [MFCC, DCT] = ma_mfcc(wavFile,p);
        [numCoeffs numFrames] = size(MFCC);        %each frame is a column-vector
        
        totalFrames = totalFrames + numFrames;
        % -------------------------------------------------------
        % Compute the FIRST and SECOND Derivatives of MFCC: Compute DYN_MFCC
        
        fd_mfcc = zeros(size(MFCC));    % First_Derivative mfcc
        sd_mfcc = zeros(size(MFCC));    % Second_Derivative mfcc
        
        % extract dyn_mfcc:: Check the formulae!!
        for frame = 1:numFrames
            if frame > 1 && frame < numFrames
                fd_mfcc(:,frame) = ((-1)^2 + (1)^2)^(-1) * ((-1)*fd_mfcc(:,frame-1)+(1)*fd_mfcc(:,frame+1));
                                sd_mfcc(:,frame) = ((1)^2 +(-2)^2 + (1)^2)^(-1) * ((1)*fd_mfcc(:,frame-1)+ (-2)*fd_mfcc(:,frame) +(1)*fd_mfcc(:,frame+1));
            elseif frame == 1
                fd_mfcc(:,frame) =  ((-1)^2 + (1)^2)^(-1)*((-1)*fd_mfcc(:,frame)+ (1)*fd_mfcc(:,frame+1));
                                sd_mfcc(:,frame) = ((1)^2 +(-2)^2 + (1)^2)^(-1) *  (0 + (-1)*fd_mfcc(:,frame) +(1)*fd_mfcc(:,frame+1));
            else
                fd_mfcc(:,frame) = ((-1)^2 + (1)^2)^(-1)*((-1)*fd_mfcc(:,frame-1) + (1)*fd_mfcc(:,frame));
                                sd_mfcc(:,frame) = ((1)^2 +(-2)^2 + (1)^2)^(-1) * ((1)*fd_mfcc(:,frame-1)+ (-1)*fd_mfcc(:,frame) + 0);
            end
        end
        
        DYN_MFCC = [MFCC; fd_mfcc; sd_mfcc];
        FD_DYN_MFCC = [MFCC; fd_mfcc];
        warning('dyn_mfcc_allSongs')
        
        
        
        % -------------------------------------------------------
        % storing MFCC and DYN_MFCC into same file.
        save(mfcc_song_fileName,'MFCC','DYN_MFCC','FD_DYN_MFCC');
        
        
    end % iterate over all songs
end % iterate over all Files

% -------------------------------------------------------------------------
% Computing Aggregate Information
[dims_FD_DYN_MFCC, frameSong] = size(DYN_MFCC);
[dims_DYN_MFCC, frameSong] = size(DYN_MFCC);
[dims_MFCC, frameSong] = size(MFCC);

aggData_mFile.TOTAL_FRAMES = totalFrames;
aggData_mFile.TOTAL_SONGS = countSong;

aggData_mFile.dims_FD_DYN_MFCC = dims_FD_DYN_MFCC;
aggData_mFile.dims_DYN_MFCC = dims_DYN_MFCC;
aggData_mFile.dims_MFCC = dims_MFCC;

warning('dyn_mfcc : has been extracted and saved for every song in DB.');


% -------------------------------------------------------------------------
% Clear the large data-variables in memory
clear all; close all;




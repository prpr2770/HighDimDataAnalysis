%{
Script to do the following:
1. Read all the songs and extract all their mfcc information. Also, obtain
vector representing genreIDs of songs. 
1.a) compute the dynamic_mfcc: first and second derivatives
2. Normalize the dyn_mfcc vectors per dimension
3. Store/archive mfcc and dyn_mfcc for all songs!

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
status = mkdir(matDataDirName)

% data directory to store mfcc and dyn_mfcc of all song
aggregateDataDirName = fullfile(tracksDirName,'dynMfcc_data_aggregate')
status = mkdir(aggregateDataDirName)

mfcc_all_songs_fileName = strcat(aggregateDataDirName,'\dyn_mfcc_allSongs.mat');
mFile = matfile(mfcc_all_songs_fileName,'Writable',true);

% read all songs and extract mfcc
Files=dir(tracksDirName);
length(Files)


% -------------------------------------------------------------------------
% p = struct('fs',11025,'visu',0,'fft_size',256,'hopsize',128,'num_ceps_coeffs',0,'dB_max',96,'mel_filt_bank','auditory-toolbox');
p = struct('fs',11025,'visu',0,'fft_size',256,'hopsize',128,'num_ceps_coeffs',20,'dB_max',96);

% mfcc_allSongs = [];
dyn_mfcc_allSongs = [];

countSong = 0;

for k=1:length(Files)
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
        [numCoeffs numFrames] = size(MFCC)         %each frame is a column-vector
        
%         % archive the MFCCs of all songs
%         mfcc_allSongs = [mfcc_allSongs MFCC];
%         warning('mfcc_allSongs')
%         size(mfcc_allSongs)

        % -------------------------------------------------------
        fd_mfcc = zeros(size(MFCC));    % First_Derivative mfcc
        sd_mfcc = zeros(size(MFCC));    % Second_Derivative mfcc
        
        % extract dyn_mfcc
        for frame = 1:numFrames
            if frame > 1 && frame < numFrames
                fd_mfcc(:,frame) = ((-1)^2 + (1)^2)^(-1) * ((-1)*fd_mfcc(:,frame-1)+(1)*fd_mfcc(:,frame+1));
                sd_mfcc(:,frame) = ((1)^2 +(-2)^2 + (1)^2)^(-1) * ((1)*fd_mfcc(:,frame-1)+ (-2)*fd_mfcc(:,frame) +(1)*fd_mfcc(:,frame+1));
            elseif frame == 1
                fd_mfcc(:,frame) =  (0 + (1)*fd_mfcc(:,frame+1));
                sd_mfcc(:,frame) = ((-1)^2+(1)^2)^(-1) * (0 + (-1)*fd_mfcc(:,frame) +(1)*fd_mfcc(:,frame+1));
            else
                fd_mfcc(:,frame) = ((1)*fd_mfcc(:,frame-1) + 0);
                sd_mfcc(:,frame) = ((1)^2  + (-1)^2)^(-1) * ((1)*fd_mfcc(:,frame-1)+ (-1)*fd_mfcc(:,frame) + 0);
            end
        end
        
        DYN_MFCC = [MFCC; fd_mfcc; sd_mfcc];
        warning('dyn_mfcc_allSongs')
        size(DYN_MFCC)
        
        save(mfcc_song_fileName,'MFCC','DYN_MFCC');
        
        % archive the dyn_mfcc_ofAllSongs!
        if countSong ~= 1
            [nrows ncols] = size(mFile,'dyn_mfcc_allSongs');
            mFile.dyn_mfcc_allSongs(:,ncols+1:ncols+numFrames) = DYN_MFCC;
        else
            mFile.dyn_mfcc_allSongs = DYN_MFCC;
        end
    
%         warning('dyn_mfcc_allSongs')
%         size(dyn_mfcc_allSongs)
    
        % -------------------------------------------------------
        % Store the mfcc and dyn_mfcc of each song separately!
        
        
    end % iterate over all songs
    
        
end % iterate over all Files

% save(mfcc_all_songs_fileName,'dyn_mfcc_allSongs');


% ----------------------------------------------------
% modify script to compute the mean, std_dev and the normalized vectors, directly
% from the .mat file. without loading the whole data into memory.

% -------------------------------------------------------
% Normalize the dyn_mfcc_allSongs:
totalFramesInBag = size(dyn_mfcc_allSongs,2);

mean_dyn_mfcc = 1/totalFramesInBag * sum(dyn_mfcc_allSongs,2); %rowsum
mean_dyn_mfcc_mat = repmat(mean_dyn_mfcc,1,totalFramesInBag);

std_dyn_mfcc = (1/totalFramesInBag * sum((dyn_mfcc_allSongs - mean_dyn_mfcc_mat).^2,2 )).^(1/2);
std_dyn_mfcc_mat = repmat(std_dyn_mfcc,1,totalFramesInBag);

dyn_mfcc_allSongs = (dyn_mfcc_allSongs - mean_dyn_mfcc_mat)./std_dyn_mfcc_mat;

% -------------------------------------------------------
% Archive mfcc_allSongs, dyn_mfcc_allSongs 

save(mfcc_all_songs_fileName, 'mean_dyn_mfcc', 'std_dyn_mfcc', 'countSong','-append');

% Archive: norm_dyn_mfcc_allSongs mean_dyn_mfcc std_dyn_mfcc


% -------------------------------------------------------
% Clear the large data-variables in memory
clear mfcc_allSongs dyn_mfcc_allSongs MFCC fd_mfcc sd_mfcc DYN_MFCC
warning('cleared: mfcc_allSongs dyn_mfcc_allSongs MFCC fd_mfcc sd_mfcc DYN_MFCC')
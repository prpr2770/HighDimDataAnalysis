% % % % read properties of the wav file
close all; clear all; clc
% create structure for using ma_mfcc
% p = struct('fs',{},'visu',0,'fft_size',256,'hopsize',128,'num_ceps_coeffs',20,'use_first_coeff',1,'mel_filt_bank','auditory-toolbox','dB_max',96);
p = struct('fs',11025,'visu',0,'fft_size',256,'hopsize',128,'num_ceps_coeffs',20,'dB_max',96,'mel_filt_bank','auditory-toolbox');

% ================================================================
% Access the truth File with information about Genres

% determine genreName from truth
dirName = 'H:\HighDimData\Project\ecen5322\Volumes\project'
truthFileName = 'ground_truth.csv'
fullFileName = fullfile(dirName,truthFileName)
fileID = fopen(fullFileName);
C = textscan(fileID,'%s %s','Delimiter',',')
fclose(fileID)

% remove quotation marks
C{1} = strrep(C{1},'"','');
C{2} = strrep(C{2},'"','');
% remove 'tracks/' prefix from the filename
C{1} = strrep(C{1},'tracks/','');
% rename '.mp3' as '.wav' in suffix
C{1} = strrep(C{1},'.mp3','.wav');
% celldisp(C)


wavfileGenreDict = containers.Map(C{1},C{2});

genreKeys = {'rock_pop', 'classical', 'electronic', 'jazz_blues', 'metal_punk', 'world'};
genreValues = {'rockp', 'class', 'elect', 'jazzb', 'metal', 'world'};
genreRenameDict = containers.Map(genreKeys, genreValues);

% ================================================================
% Accessing the tracks and updating their names

dirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\'
Files=dir(dirName)
length(Files)

% ===================================================================
%{

% -----------------------------------------------------------
% determine Number of Valid SongFiles
fileSizes = zeros(size(Files));

for i = 1:length(Files)
    fileSizes(i) = Files(i).bytes;
end

totalSongs = sum(fileSizes > 0);
totalSamplesInEachSong = zeros(1, totalSongs);
% -----------------------------------------------------------

% 1. Determine length of all songs
countSong = 0;
for k=1:length(Files)
    if (Files(k).bytes >0)
        countSong = countSong + 1;
        wavName=Files(k).name;
       
        %remove the .wav extension
        fileName = wavName(1:end-4);
        % construct full file-name
        wavFileName = fullfile(dirName,wavName);
        
        % start reading the files
        try
            % read audio-file
            wavInfo = audioinfo(wavFileName);

            % determine legnth of song
            totalSamplesInEachSong(countSong) = wavInfo.TotalSamples;
        catch E1
            msg = strcat('Error in audioread: ',wavFileName);
            warning(msg)
        end

    end
end

% -----------------------------------------------------------
% 2. Determine minimum Length of all available songs

[counts,centers] = hist(totalSamplesInEachSong);
fig1 = figure(1)
plot(centers,counts)
shg
minSongLength = min(totalSamplesInEachSong)

%}

% minSize = 2^17;
minSongLength = 2^17;
% -----------------------------------------------------------
% =======================================================================

% define new Directory Name
csvDataDir = fullfile(dirName,'csvData_cut')
status = mkdir(csvDataDir)

matDataDir = fullfile(dirName,'matData_cut')
status = mkdir(matDataDir)


fid_allSongNames = fopen('allSongNames.txt', 'at');

% =====================================================================
% =====================================================================

mfcc_all_fileName = strcat(matDataDir,'\allSongs_vec_mfcc.mat');
mfcc_all_mat_fileName = strcat(matDataDir,'\allSongs_mat_mfcc.mat');

mfcc_genre_fileName = strcat(matDataDir,'\allGenre_Songs_vec_mfcc.mat');
mfcc_genre_mat_fileName = strcat(matDataDir,'\allGenre_Songs_mat_mfcc.mat');


% =====================================================================
% =====================================================================
% Read and create MFCC Files
countSong = 0;

countGenreSongs = zeros(1, length(genreKeys));

% create txt file containing all valid song names



for k=1:length(Files)
    if (Files(k).bytes >0)
        countSong = countSong + 1
        wavName=Files(k).name;
        
        % determine genre
        genreName = wavfileGenreDict(wavName);
        renameGenre = genreRenameDict(genreName);
        
        countGenreSongs = updateGenreCount(genreKeys, countGenreSongs, genreName);
        
        
        %remove the .wav extension
        fileName = wavName(1:end-4);
        % construct full file-name
        wavFileName = fullfile(dirName,wavName);
        
        
        
        
        % start reading the files
        try
            % read audio-file
            [wavFile, Fs] = audioread(wavFileName);
            wavInfo = audioinfo(wavFileName);
            
            % -------------------------------------------------------
        catch E1
            msg = strcat('Error in audioread: ',wavFile);
            warning(msg)
        end
        
        % ============================================================
        % -------------------------------------------------------
        % truncate the sound file to extract only a segment of the
        % song - of minLength of all songs in Directory
        
        windowSize = 2^(floor(log2(minSongLength))); % can be changed later!
        halfWindow = windowSize/2;
        
        %find the central data-point:
        lenWave = length(wavFile);
        if mod(lenWave,2) ~= 0 % Odd length
            cntr = (lenWave + 1)/2;
        else
            cntr = lenWave/2;
        end
        
        
        % split the signal at the center
        s_half = wavFile(cntr+1:end);
        f_half = wavFile(1:cntr);
        
        % split the signal at the center
        s_half_window = s_half(1:halfWindow);
        f_half_window = f_half(end-halfWindow+1:end);
        
        
        % generate truncated songFile
        truncatedFile = [f_half_window; s_half_window];
        
        size(truncatedFile);
        % ============================================================
        % generate mfcc and dct
        [MFCC, DCT] = ma_mfcc(truncatedFile,p);
        
        %{
        
%         warning('Sizes of the MFCC');
       size(MFCC);
       size(DCT);
        % ================================================================
        % write CSV files
        
        extName = '_full';
        dotCsv = '.csv';
        %         mpcc_csvName = strcat(fileName,extName,num2str(1),'_mpcc',dotCsv);
        mpcc_csvName = strcat(renameGenre,'_',fileName,extName,'_mfcc',dotCsv);
        dct_csvName = strcat(renameGenre,'_',fileName,extName,'_dct',dotCsv);
        
        dataFile = fullfile(csvDataDir,mpcc_csvName);
        csvwrite(dataFile,MFCC);
        % dlmwrite(dataFile,MFCC,'delimiter',',','precision',5)
        clear dataFile;
        
        dataFile = fullfile(csvDataDir,dct_csvName);
        csvwrite(dataFile,DCT);
        clear dataFile;
        
        
        % ================================================================
        % write MAT files: Individual data files

        extName = '_full';
        %         mpcc_csvName = strcat(fileName,extName,num2str(1),'_mpcc',dotCsv);
        mpcc_matName = strcat(matDataDir,'\',renameGenre,'_',fileName,extName,'_mfcc.mat');
        dct_matName = strcat(matDataDir,'\',renameGenre,'_',fileName,extName,'_dct.mat');
        
        %         dataFile = fullfile(matDataDir,mpcc_matName);
        %         save mpcc_matName MFCC;
        save(mpcc_matName, 'MFCC');
        clear dataFile;
        
        %         dataFile = fullfile(matDataDir,dct_matName);
        %         save dct_matName DCT ;
        save(dct_matName, 'DCT');
        clear dataFile;
        %}

        
%{
        % ================================================================
        % write MAT files: Appended MFCC Matrices with first col discarded
        
        mfcc_alt = MFCC(:,2:end);  % Coeffs x Frames
        numFrames = size(mfcc_alt,2);
        genreIndex = getGenreIndex(genreKeys, genreName);
        genreIndexVec = genreIndex*ones(1,numFrames);

        % ----------------------------------------------------------------
        % Store All Songs MFCC : As MFCC Matrix
        % X: Data Frames with MFCC coefficients
        % Y: Data Labels per Frame

        if (countSong > 1)
            [nrows,ncols]=size(mfccDataFile,'X');
            numFrames = size(mfcc_alt,2);
            mfccDataFile.X(:,ncols+1:ncols+numFrames) = mfcc_alt;
            [nrows,ncols]=size(mfccDataFile,'Y');
            mfccDataFile.Y(nrows,ncols+1:ncols+numFrames) = genreIndexVec;
        else
            %             mfccDataFile = matfile(mfcc_all_fileName, 'Writable', isWritable);
            mfccDataFile = matfile(mfcc_all_mat_fileName, 'Writable', true);
            mfccDataFile.X = mfcc_alt;
            mfccDataFile.Y = genreIndexVec;        % Store Label of each song
        end

%}
        
        % ================================================================
        % Store Genre-based MFCC files
        
        mfcc_alt = MFCC(:,2:end);  % Coeffs x Frames; Discard col-1
        numFrames = size(mfcc_alt,2);
        genreIndex = getGenreIndex(genreKeys, genreName);
        
        if (sum(countGenreSongs) > 1)
            
            if(countGenreSongs(genreIndex)> 1)
                % obtain size of dataset stored for given genre
                eval(['[nrows,ncols]=size(mfccGenreFile.Genre' num2str(genreIndex) ');' ])
                % store data value
                eval(['mfccGenreFile.Genre' num2str(genreIndex) '(:,ncols+1:ncols+numFrames) = mfcc_alt;'])
            else
                eval(['mfccGenreFile.Genre' num2str(genreIndex) '= mfcc_alt;'] );
            end
            
        else
            % INITIALIZE THE GENRE DATA FILE
            %             warning('Initialize: Genre Data')
            mfccGenreFile = matfile(mfcc_genre_mat_fileName, 'Writable', true);

            % get index; update each Genre type data as Genre1, Genre2, etc.
            eval(['mfccGenreFile.Genre' num2str(genreIndex) '= mfcc_alt;'] );
        end



%{
        % ================================================================
        % write MAT files: Appended MFCC vectors
        
        mfcc_vec = reshape(MFCC,[],1); % create Column Vector
        genreIndex = getGenreIndex(genreKeys, genreName);

        % ----------------------------------------------------------------
        % Store All Songs MFCC : As MFCC Vector
        % X: Data vectors
        % Y: Data Labels

        if (countSong > 1)
            [nrows,ncols]=size(mfccDataFile,'X');
            mfccDataFile.X(:,ncols+1) = mfcc_vec;
            [nrows,ncols]=size(mfccDataFile,'Y');
            mfccDataFile.Y(nrows,ncols+1) = genreIndex;
        else
            %             mfccDataFile = matfile(mfcc_all_fileName, 'Writable', isWritable);
            mfccDataFile = matfile(mfcc_all_fileName, 'Writable', true);
            mfccDataFile.X = mfcc_vec;
            mfccDataFile.Y = genreIndex;        % Store Label of each song
        end
        % delete nrows ncols
%}
      

%{
        % ----------------------------------------------------------------
        % Store Genre-based Dataset
        
        genreIndex = getGenreIndex(genreKeys, genreName);
        %         countGenreSongs
        
        if (sum(countGenreSongs) > 1)
            
            if(countGenreSongs(genreIndex)> 1)
                %             warning('Updating with NEW SONG data')
                % obtain size of dataset stored for given genre
                eval(['[nrows,ncols]=size(mfccGenreFile.Genre' num2str(genreIndex) ');' ])
                % store data value
                eval(['mfccGenreFile.Genre' num2str(genreIndex) '(:,ncols+1) = mfcc_vec;'])
            else
                eval(['mfccGenreFile.Genre' num2str(genreIndex) '= mfcc_vec;'] );
            end
            
        else
            % INITIALIZE THE GENRE DATA FILE
            %             warning('Initialize: Genre Data')
            mfccGenreFile = matfile(mfcc_genre_fileName, 'Writable', true);
            
            %             % initialize all genreTypes Data
            %             for i=1:length(genreKeys)
            %                 eval(['mfccGenreFile.Genre' num2str(i) '= [];'] );
            %             end
            
            % get index; update each Genre type data as Genre1, Genre2, etc.
            eval(['mfccGenreFile.Genre' num2str(genreIndex) '= mfcc_vec;'] );
        end
        
        % ================================================================
        
%}
        

    end
end






fclose(fid_allSongNames);        























% % % % % % =====================================================
% % % % % % specify the wavfile destination
% % % % % hfile1 = 'artist_1_album_1_track_1.wav'
% % % % % [wavFile, Fs] = audioread(hfile1);
% % % % % hfile2 = 'artist_78_album_1_track_3.wav'
% % % % %
% % % % % try
% % % % % [wavFile, Fs] = audioread(hfile2);
% % % % % catch ME
% % % % %     switch ME.identifier
% % % % %         case 'MATLAB:audiovideo:audioread:FileTypeNotSupported'
% % % % % %              msg = ['Error in Audio file: File Type Not Supported: ', hfile2];
% % % % % %              causeException = MException('MATLAB:readMusicFile:fileAccessError',msg);
% % % % % %             ME = addCause(ME,causeException);
% % % % % %             rethrow(ME)
% % % % % %         otherwise
% % % % %             warning('Audio File Not Read!')
% % % % %     end
% % % % % end
% %
% % % % %
% % % % % % gather other info regarding the audio-file
% % % % % wavInfo = audioinfo(hfile1)
% % % % %
% % % % % % determine time-duration of the file.
% % % % % trueDuration = wavInfo.Duration
% % % % % duration = numel(wavFile) / Fs % gives the number of seconds of data
% % % % %
% % % % % % determine mfcc
% % % % % [mfcc1, DCT] = ma_mfcc(wavFile,p);
% % % % %
% % % % % % =====================================================
% % % % % % specify the wavfile destination
% % % % % hfile2 = 'artist_78_album_1_track_3.wav'
% % % % % [wavFile, Fs] = audioread(hfile2);
% % % % % wavInfo = audioinfo(hfile2)
% % % % % [mfcc2, DCT] = ma_mfcc(wavFile,p);
% % % % %
% % % % % % =====================================================
% % % % % % determine if the mfcc are equal or not?
% % % % % isequal(mfcc1,mfcc2)
% % % % %
% % % % % % =====================================================
% % % % % % store mpcc and dct in a format accessible from python
% % % % % % Format: CSV
% %
% %
% %

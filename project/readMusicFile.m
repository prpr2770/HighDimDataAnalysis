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


wavfileGenreDict = containers.Map(C{1},C{2})

genreKeys = {'rock_pop', 'classical', 'electronic', 'jazz_blues', 'metal_punk', 'world'};
genreValues = {'rockp', 'class', 'elect', 'jazzb', 'metal', 'world'};
genreRenameDict = containers.Map(genreKeys, genreValues);

% ================================================================
% Accessing the tracks and updating their names

dirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\'
Files=dir(dirName)

% define new Directory Name
csvDataDir = fullfile(dirName,'csvData');
status = mkdir(csvDataDir)

for k=1:length(Files)
    if (Files(k).bytes >0)
        wavName=Files(k).name;
        
        % determine genre
        genreName = wavfileGenreDict(wavName)
        renameGenre = genreRenameDict(genreName)
        
        %remove the .wav extension
        fileName = wavName(1:end-4);
        % construct full file-name
        wavFileName = fullfile(dirName,wavName);
        
        % start reading the files
        try
            % read audio-file
            [wavFile, Fs] = audioread(wavFileName);
            wavInfo = audioinfo(wavFileName);
            
            % generate mfcc and dct
            [MFCC, DCT] = ma_mfcc(wavFile,p);
        catch E1
            msg = strcat('Error in audioread: ',wavFile);
             warning(msg)
        end
            

        extName = '_full';
        dotCsv = '.csv';
%         mpcc_csvName = strcat(fileName,extName,num2str(1),'_mpcc',dotCsv);
        mpcc_csvName = strcat(renameGenre,'_',fileName,extName,'_mpcc',dotCsv);
        dct_csvName = strcat(renameGenre,'_',fileName,extName,'_dct',dotCsv);
        
        dataFile = fullfile(csvDataDir,mpcc_csvName);
        csvwrite(dataFile,MFCC);
        % dlmwrite(dataFile,MFCC,'delimiter',',','precision',5)
        clear dataFile
        
        dataFile = fullfile(csvDataDir,dct_csvName);
        csvwrite(dataFile,DCT);
        clear dataFile

    end
end






% % % % =====================================================
% % % % specify the wavfile destination
% % % hfile1 = 'artist_1_album_1_track_1.wav'
% % % [wavFile, Fs] = audioread(hfile1);
% % % hfile2 = 'artist_78_album_1_track_3.wav'
% % % 
% % % try
% % % [wavFile, Fs] = audioread(hfile2);
% % % catch ME
% % %     switch ME.identifier
% % %         case 'MATLAB:audiovideo:audioread:FileTypeNotSupported'
% % % %              msg = ['Error in Audio file: File Type Not Supported: ', hfile2];
% % % %              causeException = MException('MATLAB:readMusicFile:fileAccessError',msg);
% % % %             ME = addCause(ME,causeException);
% % % %             rethrow(ME)
% % % %         otherwise
% % %             warning('Audio File Not Read!')
% % %     end
% % % end

% % % 
% % % % gather other info regarding the audio-file
% % % wavInfo = audioinfo(hfile1)
% % % 
% % % % determine time-duration of the file.
% % % trueDuration = wavInfo.Duration
% % % duration = numel(wavFile) / Fs % gives the number of seconds of data
% % % 
% % % % determine mfcc
% % % [mfcc1, DCT] = ma_mfcc(wavFile,p);
% % % 
% % % % =====================================================
% % % % specify the wavfile destination
% % % hfile2 = 'artist_78_album_1_track_3.wav'
% % % [wavFile, Fs] = audioread(hfile2);
% % % wavInfo = audioinfo(hfile2)
% % % [mfcc2, DCT] = ma_mfcc(wavFile,p);
% % % 
% % % % =====================================================
% % % % determine if the mfcc are equal or not?
% % % isequal(mfcc1,mfcc2) 
% % % 
% % % % =====================================================
% % % % store mpcc and dct in a format accessible from python
% % % % Format: CSV




% Code script: Reads the audio-files from directory, and creates windows/snippets

% % % % read properties of the wav file
close all; clear all; clc

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

dirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
Files=dir(dirName);

% define new Directory Name
csvDataDir = fullfile(dirName,'csvSnippetsData')
status = mkdir(csvDataDir)

for k=1:length(Files)
    k
    if (Files(k).bytes >0)
        wavName=Files(k).name;
        
        % determine genre
        genreName = wavfileGenreDict(wavName);
        renameGenre = genreRenameDict(genreName);
        
        %------------------------------------------------------------------
        % Read the Audio Files
        
        %remove the .wav extension
        fileName = wavName(1:end-4);
        % construct full file-name
        wavFileName = fullfile(dirName,wavName);
        
        % start reading the files
        try
            % read audio-file
            [wavFile, Fs] = audioread(wavFileName);
            wavInfo = audioinfo(wavFileName);
        catch E1
            msg = strcat('Error in audioread: ',wavFile);
             warning(msg)
        end
        %-----------------------------------------------------------------
        % Extract windows from the song
        
        % ensure the song is a row-vector
        if size(wavFile,1) ~= 1
            wavFile = wavFile';
        end
        
        wavWindows = extractWindowsFromSong(wavFile, Fs);
        
        % archive/ save the snippets as .mat file
        snippet_dataName = strcat(renameGenre,'_',fileName,'_windows','.csv');
        dataFile = fullfile(csvDataDir,snippet_dataName);
        csvwrite(dataFile, wavWindows);
        clear dataFile
      
        %-----------------------------------------------------------------        
        % generate mfcc and dct for the song-snippets

        % Determine number of windows generated
        numWindows = size(wavWindows,1);
        
        for win = 1:numWindows
            % obtain snippet
            snippet = wavWindows(win,:);
            
            % create structure for using ma_mfcc
            % p = struct('fs',{},'visu',0,'fft_size',256,'hopsize',128,'num_ceps_coeffs',20,'use_first_coeff',1,'mel_filt_bank','auditory-toolbox','dB_max',96);
            p = struct('fs',11025,'visu',0,'fft_size',512,'hopsize',128,'num_ceps_coeffs',20,'dB_max',96,'mel_filt_bank','auditory-toolbox');
            
            % gnerate mfcc, dct
            [MFCC, DCT] = ma_mfcc(snippet,p);
            
            extName = strcat('_window_',num2str(win));
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
        end % end iterating over windows for given song.

    end
end


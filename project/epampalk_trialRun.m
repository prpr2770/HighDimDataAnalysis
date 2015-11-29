%{

Script that does the following
1. Reads track filenames and enumerates it into a .txt file
2. songFileNames.txt is input into FeatureExtraction function
3. Distance_Matrix is computed for all the Songs. 

4. Determine Genres of the songs from truthFile and save a .mat file 
containing genreIndex for the songs.
%}

clear all;close all; clc

% ================================================================
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

% genreKeys = {'rock_pop', 'classical', 'electronic', 'jazz_blues', 'metal_punk', 'world'};
genreKeys = unique(wavfileGenreDict.values())
% create Dictionary of GenreName with IndexKey
numGenres = length(genreKeys);
genreIndx = 1:numGenres;
genreIndxDict = containers.Map(genreKeys, genreIndx);
% ================================================================


tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\'
Files=dir(tracksDirName)
length(Files)

% define new Directory Name - to store Features
g1c_SongFeatures_Dir = fullfile(tracksDirName ,'g1c_SongFeatures\')
status = mkdir(g1c_SongFeatures_Dir)


% -------------------------------------------
% create file to store the songNames
allSongsListName = 'allSongNames.txt';
allSongsListFileName = fullfile('H:\HighDimData\Project\ecen5322\Volumes\project\',allSongsListName);
% if file exists delete it!
if (exist(allSongsListFileName, 'file') == 2)
    delete(allSongsListFileName);
end
% open File
fid_allSongNames = fopen(allSongsListFileName, 'at');



% -------------------------------------------
% track total valid songs
countSong = 0;

for k=1:length(Files)
    if (Files(k).bytes >0)
        countSong = countSong + 1;
        wavName=Files(k).name;
        
        % determine genre
        genreName = wavfileGenreDict(wavName);
        genreIndex = genreIndxDict(genreName);
        Files(k).genre = genreName;
        Files(k).genreIndex = genreIndex;
        
        % determine full filename
        fullFileName = fullfile(tracksDirName,wavName);
        % append into text file
        if countSong > 1
            fprintf(fid_allSongNames, '\n%s', fullFileName);
        else 
            countSong
            disp(fullFileName)
            fprintf(fid_allSongNames, '%s', fullFileName);
            
        end
        
    end
end

fclose(fid_allSongNames)

songGenres = zeros(1,countSong);
countSong = 0;
for k=1:length(Files)
    if (Files(k).bytes >0)
        countSong = countSong + 1;
        songGenres(countSong) = Files(k).genreIndex;

    end
end

songGenres_filename = 'songGenres.mat';
songGenres_full_filename = fullfile(tracksDirName,'g1c_SongFeatures\', songGenres_filename) ;
save(songGenres_full_filename,'songGenres','genreKeys','wavfileGenreDict');


%{
% ==========================================================================
% allSongsListFileName : contains list of all songs to be processed.

% function FeatureExtraction(in_file,out_dir): MA_Toolbox (Elias Pampalk)
FeatureExtraction(allSongsListFileName,g1c_SongFeatures_Dir);

% ==========================================================================
% Compute Distance Matrix
distMatrix_filename = 'distance_matrix.txt';
distMatrix_full_filename = fullfile(tracksDirName,'g1c_SongFeatures\', distMatrix_filename) ;
% ComputeSimilarities(out_dir,'somepath/distance_matrix.txt')
ComputeSimilarities(g1c_SongFeatures_Dir,distMatrix_full_filename);
%}





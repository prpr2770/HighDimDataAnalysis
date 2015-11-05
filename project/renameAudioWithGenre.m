% Matlab code to rename all the files in a particular directory to
% incorporate the GENRE they belong to, which is given in another file.

close all; clear all; clc

% ================================================================
% Access the truth File

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

% take one file
k = 10
wavName=Files(k).name % 'song1_track2_part1.wav'

% determine genre
genreName = wavfileGenreDict(wavName)
renameGenre = genreRenameDict(genreName)


% % 
% % % rename with genre prefix
% % newWavName = strcat(genreName,wavName);
% % 
% % % how to rename a file in Matlab
% % movefile('oldname.m','newname.m')
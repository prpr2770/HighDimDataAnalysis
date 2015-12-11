%{

Script to do the following:
+ read genreKeys
+ read genreTestSongs file

create genreSongs vector;
%}

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\testTracks\';

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'aggregate_song_data')


% -------------------------------------------------------
% file 1: 
genreSongs_fileName = strcat(aggregateDataDirName,'\testTracks_genres.txt');

formatSpec = '%s \t %s'
fileID = fopen(genreSongs_fileName,'r')
C = textscan(fileID,formatSpec, 'delimiter', '\t');
fclose(fileID);
% celldisp(C);

% -----------------------------------------------------
% obtain the songGenre.mat

% matfile to store aggData_allSongs
genreKeys_fileName = strcat(aggregateDataDirName,'\songGenres.mat');
genreKeys_mFile = matfile(genreKeys_fileName,'Writable', true);

GENRE_KEYS = genreKeys_mFile.genreKeys;

totalSongs = length(C{2});

TESTSONG_GENRES = zeros(1,totalSongs);

for i = 1:totalSongs
    songGenre = C{2}{i};
    
    switch songGenre
        case 'classical'
            TESTSONG_GENRES(i) = 1;
            
        case 'electronic'
            TESTSONG_GENRES(i) = 2;
        case 'jazz'
            TESTSONG_GENRES(i) = 3;
            
        case 'punk'
            TESTSONG_GENRES(i) = 4;
        case 'rock'
            TESTSONG_GENRES(i) = 5;
        case 'world'
            TESTSONG_GENRES(i) = 6;
           
    end
    
end

genreKeys_mFile.TESTSONG_GENRES = TESTSONG_GENRES;

% BAM_make_confusion_matrix :: code to generate the confusion Matrix.
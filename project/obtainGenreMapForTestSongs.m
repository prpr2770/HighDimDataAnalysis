tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\testTracksFolders\';

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'aggregate_song_data')

genreSongs_fileName_csv = strcat('H:\HighDimData\Project\ecen5322\Volumes\project','\testTracks_genres.txt');
fid=fopen(genreSongs_fileName_csv,'wt');

Files = dir(tracksDirName)

songCount = 0;
genreKeys = [];
genreCount = 0;

TEST_GENRE_ID = [];
for file = 1:length(Files)
    if (~strcmp(Files(file).name,'.') && ~strcmp(Files(file).name,'..'))
        folderName = Files(file).name;
        newTracksDirName = strcat(tracksDirName,folderName,'\');
        Tracks = dir(newTracksDirName);
        
        genreCount = genreCount + 1;
        genreKeys = [genreKeys folderName];
        
        for song = 1: length(Tracks)
            if (Tracks(song).bytes >0 )
                songCount = songCount + 1;
                wavfileName=Tracks(song).name;
                
                
                msg = sprintf('%s \t %s\n',wavfileName,folderName);
                fprintf(fid,msg);
                
             
                
            end
            
            
            
        end
        
        
    end
    
end



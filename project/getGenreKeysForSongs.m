function [genreKeys songGenres] = getGenreKeysForSongs(tracksDirName)
% Determine Genres for each type: read the genreFile to determine the index 
songGenres_filename = 'songGenres.mat';
songGenres_full_filename = fullfile(tracksDirName,'g1c_SongFeatures\', songGenres_filename) ;
% save(songGenres_full_filename,'songGenres','genreKeys','wavfileGenreDict');
songGenresData = load(songGenres_full_filename);

songGenres = songGenresData.songGenres;
genreKeys = songGenresData.genreKeys;
end


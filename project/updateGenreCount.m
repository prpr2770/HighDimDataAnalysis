function countSongs =  updateGenreCount(genreKeys, countGenreSongs, genreName)
%{
function to update genreCount: number of songs of each Genre detected in DB
genreKeys = {'rock_pop', 'classical', 'electronic', 'jazz_blues', 'metal_punk', 'world'};

countGenreSongs : vector with the counts
genreName : genre of current song.
%}

indx = getGenreIndex(genreKeys, genreName);

numGenres = length(genreKeys);

% update the genreIndicator vector
genreIndicator = zeros(1,numGenres);
genreIndicator(indx) = 1;

% update genreCounts
countSongs = countGenreSongs + genreIndicator;

end
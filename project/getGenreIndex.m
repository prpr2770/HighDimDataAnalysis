function indx = getGenreIndex(genreKeys, genreName)

% create a dictionary of Genres and their Indices
numGenres = length(genreKeys);
genreIndx = 1:numGenres;
genreIndxDict = containers.Map(genreKeys, genreIndx);

% determine genreIndex for the given genre
indx = genreIndxDict(genreName);

end
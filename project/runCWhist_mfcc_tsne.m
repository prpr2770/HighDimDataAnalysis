%{
Script that does the following:

1. Read all the cw_hist for each song. 
2. Implement ITML on the code.
3. Execute tSNE and visualize. 	

%}
% =========================================================================

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'dynMfcc_data')

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'dynMfcc_data_aggregate')


% -------------------------------------------------------------------------------

% matfile: storing 'CODEWORD_HIST_ALLSONGS'
cwHist_fileName = strcat(aggregateDataDirName,'\cw_histogram_mfcc_allSongs.mat');
cwHist_mFile = matfile(cwHist_fileName);

% matfile: storing 'CODEWORD_HIST_ALLSONGS'
cwDistance_fileName = strcat(aggregateDataDirName,'\cw_distance_mfcc_allSongs.mat');
cwDistance_mFile = matfile(cwDistance_fileName);

% ===================================================================================
% Reading data from files

[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

CW_HIST_DATA = cwHist_mFile.CODEWORD_HIST_ALLSONGS;

CW_HIST_DATA = CW_HIST_DATA.^(0.5);

size(CW_HIST_DATA)
size(songGenres)

% ===================================================================================
% implement tSNE on Data

perplexity = 30;
out_dims = 3;
initial_dims = 30;

fig3 = figure(3)
zdata = tsne(CW_HIST_DATA', songGenres', out_dims, initial_dims, perplexity);


% scatter plot for all the data
fig13 = figure(13)
% colormap(jet(length(genreKeys)));
colormap(prism(length(genreKeys)));

for i = 1:length(genreKeys)
% for i = 2:length(genreKeys)-1
    songsOfAGenre = find(songGenres == i);
    hue = (length(genreKeys)+1-i)*ones(length(songsOfAGenre),1);
    rad = 30*ones(length(songsOfAGenre),1);
    h = scatter3(zdata(songsOfAGenre,1),zdata(songsOfAGenre,2),zdata(songsOfAGenre,3),rad,hue,'filled')
    hold on
end
legend(genreKeys)
hold off


%{
Script that does the following:

1. Read all the cw_hist for each song. 
2. Implement ITML on the code.
3. Execute tSNE and visualize. 	

%}


clear all; close all; clc;

% =========================================================================
clear all; close all; clc;

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';
[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'dynMfcc_data')

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'dynMfcc_data_aggregate')


% -------------------------------------------------------------------------------

% matfile: storing 'CODEWORD_HIST_ALLSONGS'
cwHist_fileName = strcat(aggregateDataDirName,'\cw_histogram_allSongs.mat');
cwHist_mFile = matfile(cwHist_fileName);

% matfile: storing 'CODEWORD_HIST_ALLSONGS'
cwDistance_fileName = strcat(aggregateDataDirName,'\cw_distance_allSongs.mat');
cwDistance_mFile = matfile(cwDistance_fileName);


% ===================================================================================
% Reading data from files

[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

CW_HIST_DATA = cwHist_mFile.CODEWORD_HIST_ALLSONGS;

size(CW_HIST_DATA)
size(songGenres)
% ===================================================================================
% Implement ITML


y = songGenres'; 						% COL-VEC: songGenres
X = CW_HIST_DATA';									% ROW-VEC: cw_hist_perSong
distance_metric = MetricLearningAutotuneKnn(@ItmlAlg, y, X);

cwDistance_mFile.DIST_METRIC = distance_metric;

imagesc(distance_metric)
colorbar

% ===================================================================================
% compute Distance_Matrix

totalSongs = length(songGenres);
D = zeros(totalSongs,totalSongs);
for i = 1:totalSongs
    for j = i:totalSongs
    x_song = X(i,:);
    y_song = X(j,:);
    d_ij = (x_song - y_song)*distance_metric*(x_song - y_song)';
    if (d_ij >= 0)
        D(i,j) = d_ij.^(1/2);
    else
        warning('distance are complex')
    end
    
    end
end
Distance_Matrix = D + D';

cwDistance_mFile.DIST_MATRIX = Distance_Matrix;

% ===================================================================================
% implement tSNE on Distance Matrix

perplexity = 30;
out_dims = 3;
fig3 = figure(3)
zdata = tsne(Distance_Matrix, songGenres, out_dims, perplexity);

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

fig4 = figure(14)
imagesc(Distance_Matrix)

fig5 = figure(15)
imagesc(distance_metric)
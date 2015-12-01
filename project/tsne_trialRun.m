% Extract the mfcc data features and implement the TSNE algorithm on the same. 

clear all;close all; clc
% dir containing tracks
tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\'
% dir with song features
g1c_SongFeatures_Dir = fullfile(tracksDirName ,'g1c_SongFeatures\')


% -------------------------------------------------------------------
% reading the distMatrix.mat file

distMatrix_mat = 'distance_matrix.mat';

% reading the .mat file
matfile_distMatrix = fullfile(g1c_SongFeatures_Dir, distMatrix_mat);
distance_matrix = load(matfile_distMatrix);

% -------------------------------------------------------------------
% read the G1C_features
features_mat = 'G1C_features.mat';

% reading the .mat file
matfile_features = fullfile(g1c_SongFeatures_Dir, features_mat);
g1c_data = load(matfile_features);

% Extract the Single Gaussian Representation of Each Song!
gaussianRept_mfcc_features = g1c_data.data.feat.g1;

% -------------------------------------------------------------------

% Extract genreIDs of the songs
[genreKeys songGenres] = getGenreKeysForSongs(tracksDirName);

size(songGenres);
totalSongs = size(songGenres,2);
% -------------------------------------------------------------------
% Prepare Test-Training Data to discriminate between Classical Songs.

% mean mfcc
mean_mfcc = gaussianRept_mfcc_features.m;
% covariance
co_mfcc = gaussianRept_mfcc_features.co;
[numSongs p q] = size(co_mfcc);
co_mfcc = reshape(co_mfcc, numSongs, []);
% inverse co-variance
ico_mfcc = gaussianRept_mfcc_features.ico;

%{
% ------------------------------------------------------
% Remove songs belonging to world and classical
class_and_world_song_indxs = [find(songGenres == 1 )  find(songGenres == 6 )];

all_other_song_indxs = setdiff(1:totalSongs, class_and_world_song_indxs);

size(class_and_world_song_indxs)
size(all_other_song_indxs)



mean_mfcc = mean_mfcc(all_other_song_indxs,:);
co_mfcc = co_mfcc(all_other_song_indxs,:);
songGenres = songGenres(all_other_song_indxs)
%}


% -------------------------------------------------------------
% create X vector
methodOfX = 'mean_and_norm_covariance_mfcc';

if size(mean_mfcc,1) == size(co_mfcc,1)
    switch methodOfX
        case 'mean_and_covariance_mfcc'
            X = [mean_mfcc co_mfcc];
        
        case 'mean_and_norm_covariance_mfcc'

            % determine normalized cov_mfcc vectors
            sum_vec = sum(co_mfcc,2); %rowsum
            [rows cols] = size(co_mfcc);
            sum_mat = repmat(sum_vec,1,cols);
            norm_co_mfcc = 0.5*co_mfcc./sum_mat;
            
            % concatenate the two matrices
            X = [mean_mfcc norm_co_mfcc];
            
        case 'normalized_mean_and_covariance_mfcc'
            
            % determine normalized mean_mfcc vectors
            sum_vec = sum(mean_mfcc,2); %rowsum
            [rows cols] = size(mean_mfcc);
            sum_mat = repmat(sum_vec,1,cols);
            norm_mean_mfcc = mean_mfcc./sum_mat;
            
            % determine normalized cov_mfcc vectors
            sum_vec = sum(co_mfcc,2); %rowsum
            [rows cols] = size(co_mfcc);
            sum_mat = repmat(sum_vec,1,cols);
            norm_co_mfcc = co_mfcc./sum_mat;
            
            % concatenate the two matrices
            X = [norm_mean_mfcc norm_co_mfcc];

        case 'normalized_mean_and_zero_covariance_mfcc'
            % determine normalized mean_mfcc vectors
            sum_vec = sum(mean_mfcc,2); %rowsum
            [rows cols] = size(mean_mfcc);
            sum_mat = repmat(sum_vec,1,cols);
            norm_mean_mfcc = mean_mfcc./sum_mat;
            
            % determine normalized cov_mfcc vectors
            sum_vec = sum(co_mfcc,2); %rowsum
            [rows cols] = size(co_mfcc);
            sum_mat = repmat(sum_vec,1,cols);
            norm_co_mfcc = 0.5*co_mfcc./sum_mat;
            
            % concatenate the two matrices
            X = [norm_mean_mfcc norm_co_mfcc];
            
            
    end
    % X = mean_mfcc;
end



%{

% -----------------------------------------------------------------------
% implement tSNE
perplexity = 30;
out_dims = 3;
initial_dims = 30;

fig1 = figure(1)
ydata = tsne(X, songGenres, out_dims, initial_dims, perplexity);

% scatter plot for all the data
fig11 = figure(11)
% colormap(jet(length(genreKeys)));
colormap(prism(length(genreKeys)));

for i = 1:length(genreKeys)
% for i = 2:length(genreKeys)-1
    songsOfAGenre = find(songGenres == i);
    hue = (length(genreKeys)+1-i)*ones(length(songsOfAGenre),1);
    rad = 30*ones(length(songsOfAGenre),1);
    h = scatter3(ydata(songsOfAGenre,1),ydata(songsOfAGenre,2),ydata(songsOfAGenre,3),rad,hue,'filled')
    hold on
end
legend(genreKeys)
hold off

% -------------------------------------------------------------------
% running tSNE on the Distance Metric 

MutualDistance = distance_matrix.D;

% implement tSNE on Distance Matrix
perplexity = 30;
out_dims = 3;
fig2 = figure(2)
zdata = tsne(MutualDistance, songGenres, out_dims, initial_dims, perplexity);

% scatter plot for all the data
fig12 = figure(12)
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


%}

% =========================================================================
% Implement Manifold Learning Toolbox: SVM to detect classical or not.
% options=ml_options('Kernel','rbf','KernelParam',0.5,'NN',6);
% classifier=ml_train(X,Y,options,'lapsvm');
% 
% [real_output, binary_labels, error_rate]=ml_test(classifier, Xtest, Ytest);


% =========================================================================
% Distance Metric Learning

% X = mean_mfcc;
y = songGenres';
distance_metric = MetricLearningAutotuneKnn(@ItmlAlg, y, X);


imagesc(distance_metric)
colorbar
% ------------------------------------------------------
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

% -----------------------------------------------
% Save the computed Distance Matrix

distMatrix_mat = 'itml_distance_matrix.mat';
distMatrix_mat = strcat(methodOfX,distMatrix_mat);

% opening the .mat file
matfile_distMatrix = fullfile(g1c_SongFeatures_Dir, distMatrix_mat);

ITML_D = Distance_Matrix;
save(matfile_distMatrix,'ITML_D','distance_metric');
% ------------------------------------------------------
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
imagesc()
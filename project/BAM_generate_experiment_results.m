

%% BAM_generate_experiment_results

% input variables...
%  - anTestData (d x num_test_samples)
%  - anTrainingData (d x num_training_samples)
%    (both above are the data vectors in whatever dimensions desired)
%
%  - anTestGenreIndex (1 x num_test_samples)
%  - anTrainingGenreIndex (1 x num_training_samples)
%    (both above are lists of genre index (1-6)
%
%  PARAMETERS...
%  - num_nbrs = the number of neighbors
%  - weight_type = (1-3) for distance weighting type
%  - weight_sigma = for Gaussian distance weighting
%  - test_equal_training_data = set to 1 if test data is the same as
%                               training data
%  - 
%
% OUTPUTS...
%  - anTestData_ConfusionMatrix
%  - anTestData_PercentCorrect
%  - anTestData_AllPercentCorrect
%
%  - 
%


[nTestDim nTestSamples] = size(anTestData);
[nTrainingDim nTrainingSample] = size(anTrainingData);


anTestGengreNumSongs = zeros(1,6);
anTrainingGengreNumSongs = zeros(1,6);

for i=1:6
    anTestGenreNumSongs(i) = sum(anTestGenreIndex == i);
    anTrainingGenreNumSongs(i) = sum(anTrainingGenreIndex == i);
end




%% create complete data set of custom feature vector 1
%% or load... here just use as anData


%% create distance matrix of songs
%% run knnsearch on these songs

num_nbrs = 7;
weight_type = 3;  % '1' = 1 / (1 + dist)
                  % '2' = 1 / (1 + dist^2)
                  % '3' = exp(-(dist^2)/sigma)
                  
weight_sigma = 600; % sigma in case of weight_type = 3 above

genre_vote_weight = 1 ./ anTrainingGenreNumSongs;

test_equal_training_data = 1;  % set to 1 if training and test data are same

if test_equal_training_data
    num_nbrs_actual = num_nbrs + 1;
else 
    num_nbrs_actual = num_nbrs;
end



anTestData_classify = zeros(1,nTestSamples);
anTestData_classify_wtdist = zeros(1,nTestSamples);
anTestData_classify_wtgenre = zeros(1,nTestSamples);
anTestData_classify_wtdistgenre = zeros(1,nTestSamples);

anTestData_nbrs = zeros(num_nbrs, nTestSamples);
anTestData_nbrs_dist = zeros(num_nbrs, nTestSamples);
anTestData_genre_prob = zeros(6, nTestSamples);
anTestData_genre_prob_wtdist = zeros(6, nTestSamples);
anTestData_genre_prob_wtgenre = zeros(6, nTestSamples);
anTestData_genre_prob_wtdistgenre = zeros(6, nTestSamples);


for i=1:nTestSamples
    
    
    [nbr_index, nbr_dist] = knnsearch(anTrainingData', anTestData(:,i)','k', num_nbrs_actual, 'distance', 'minkowski');
    
    % if test data equals training data then
    % remove the first neighbor because it should be itself
    if test_equal_training_data
        nbr_index = nbr_index(2:length(nbr_index));
        nbr_dist = nbr_dist(2:length(nbr_dist));
    end
    
    
    anTestData_nbrs(:,i) = nbr_index;
    anTestData_nbrs_dist(:,i) = nbr_dist;
    
    nbr_vote = zeros(1,length(nbr_dist)) + 1;
    nbr_vote_wt = zeros(1,length(nbr_dist)) + 1;
    
    if weight_type == 1
        nbr_vote_wt = 1 ./ (1 + nbr_dist);
    end
    
    if weight_type == 2
        nbr_vote_wt = 1 ./ (1 + nbr_dist.^2);
    end
    
    if weight_type ==3
        nbr_vote_wt = exp(-(nbr_dist.^2) ./ weight_sigma);
    end
    
    for j=1:num_nbrs
        nIndex = anTrainingGenreIndex(nbr_index(j));
        
        anTestData_genre_prob(nIndex,i) = ...
         anTestData_genre_prob(nIndex,i) + 1;
        
        anTestData_genre_prob_wtdist(nIndex,i) = ...
         anTestData_genre_prob_wtdist(nIndex,i) + nbr_vote_wt(j);
     
        anTestData_genre_prob_wtgenre(nIndex,i) = ...
         anTestData_genre_prob_wtgenre(nIndex,i) + genre_vote_weight(nIndex);
        
        anTestData_genre_prob_wtdistgenre(nIndex,i) = ...
         anTestData_genre_prob_wtdistgenre(nIndex,i) + ...
         nbr_vote_wt(j) * genre_vote_weight(nIndex);
     
    end
    
    [nMaxVal nMaxIndex] = max(anTestData_genre_prob(:,i));
    anTestData_classify(i) = nMaxIndex;

    [nMaxVal nMaxIndex] = max(anTestData_genre_prob_wtdist(:,i));
    anTestData_classify_wtdist(i) = nMaxIndex;
    
    [nMaxVal nMaxIndex] = max(anTestData_genre_prob_wtgenre(:,i));
    anTestData_classify_wtgenre(i) = nMaxIndex;

    [nMaxVal nMaxIndex] = max(anTestData_genre_prob_wtdistgenre(:,i));
    anTestData_classify_wtdistgenre(i) = nMaxIndex;

end

%% make confusion matrix
anTestData_ConfusionMatrix = zeros(6,6,4);
anTestData_PercentCorrect = zeros(1,6,4);
anTestData_AllPercentCorrect = zeros(1,4);


for i=1:nTestSamples

    nTrueGenreIndex = anTestGenreIndex(i);
    
    nClassifyGenreIndex = anTestData_classify(i);
    anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 1) = ...
     anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 1) + 1;
 
    nClassifyGenreIndex = anTestData_classify_wtdist(i);
    anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 2) = ...
     anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 2) + 1;
    
    nClassifyGenreIndex = anTestData_classify_wtgenre(i);
    anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 3) = ...
     anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 3) + 1;
    
    nClassifyGenreIndex = anTestData_classify_wtdistgenre(i);
    anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 4) = ...
     anTestData_ConfusionMatrix(nClassifyGenreIndex, nTrueGenreIndex, 4) + 1;
 
end



for i=1:4
    anTestData_PercentCorrect(1,:,i) = diag(anTestData_ConfusionMatrix(:,:,i));
    anTestData_PercentCorrect(1,:,i) = anTestData_PercentCorrect(1,:,i) ./ anTestGenreNumSongs;
    
    anTestData_AllPercentCorrect(i) = sum(diag(anTestData_ConfusionMatrix(:,:,i))) / sum(anTestGenreNumSongs);
    
    disp(' ');
    disp(anTestData_ConfusionMatrix(:,:,i));
    disp('=========================================');
    disp(anTestData_PercentCorrect(1,:,i));
    disp(anTestData_AllPercentCorrect(i));
    disp(' ');
    
end

disp(' ');
disp(anTestData_AllPercentCorrect);


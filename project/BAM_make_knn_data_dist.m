


function [anTestData_nbr_dist_wt anTestData_label_prob] = BAM_make_knn_data_dist(...
    anTestData_nbr_index, anTestData_nbr_dist, nNumLabel, anTrainingGenreIndex, ...
    DistanceWeightOn, DistanceWeightType, DistanceWeightFactor, DistanceFilterThreshold, ...
    LabelWeightOn, LabelWeightType, LabelWeightFactor)
    
    % anTestData_nbr_index
    % anTestData_nbr_dist
    
    % nNumLabel
    % anTrainingGenreIndex
    
    % anTestData_nbr_index
    % anTestData_nbr_dist
        
    % DistanceWeightOn
    % DistanceWeightType
    % DistanceWeightFactor
    % DistanceFilterThreshold - set to Inf if not desired
    
    % LabelWeightOn
    % LabelWeightType
    % LabelWeightFactor
    
    
    [nNumNbrs nTestSamples] = size(anTestData_nbr_dist);
    anTestData_nbr_dist_wt = zeros(nNumNbrs,nTestSamples);
    anTestData_label_prob = zeros(nNumLabel,nTestSamples);
    
        
    % calculating the label (genre) weighting if desired
    if LabelWeightOn

        anTrainingGengreNumSongs = zeros(1,nNumLabel);

        for i=1:nNumLabel
            anTrainingGenreNumSongs(i) = sum(anTrainingGenreIndex == i);
        end
        
        
        % default label weight, if on
        genre_vote_weight = 1 ./ anTrainingGenreNumSongs;
        
        if LabelWeightType == 1
            genre_vote_weight = 1 ./ anTrainingGenreNumSongs.^LabelWeightFactor;
        end
        
        if LabelWeightType == 2  % test LabelWeightFactor = 0.05, 0.075, 0.1, 0.5
            genre_vote_weight = 1 ./ log(anTrainingGenreNumSongs * LabelWeightFactor); 
        end
                       
        
    else
        genre_vote_weight = zeros(1,nNumLabel) + 1;
    end
    
    
    
    %weight_type = 1;  % '1' = 1 / (1 + dist)
                      % '2' = 1 / (1 + dist^2)
                      % '3' = exp(-(dist^2)/sigma)

    %DistanceWeightFactor = 600; % sigma in case of weight_type = 3 above


    for i=1:nTestSamples

        nbr_index = anTestData_nbr_index(:,i);
        nbr_dist = anTestData_nbr_dist(:,i);

        nbr_vote_wt = zeros(1,length(nbr_dist)) + 1;

        
        % calculating the distance weighting if desired
        if DistanceWeightOn
        
            if DistanceWeightType == 1
                nbr_vote_wt = 1 ./ (1 + nbr_dist);
            end

            if DistanceWeightType == 2
                nbr_vote_wt = 1 ./ (1 + nbr_dist.^2);
            end

            if DistanceWeightType ==3
                nbr_vote_wt = exp(-(nbr_dist.^2) ./ DistanceWeightFactor);
            end
            
        end
        
        
        
        
        % adding the votes;
        for j=1:nNumNbrs

            if nbr_dist(j) < DistanceFilterThreshold

                nIndex = anTrainingGenreIndex(nbr_index(j));

                nThisVote = 1 * nbr_vote_wt(j) * genre_vote_weight(nIndex);                              
                
                anTestData_label_prob(nIndex,i) = ...
                 anTestData_label_prob(nIndex,i) + nThisVote;

             end
        end
    end
end






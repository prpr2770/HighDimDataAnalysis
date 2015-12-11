

function [anTestData_nbr_index anTestData_nbr_dist] = BAM_get_knn_data_info(anTrainingData, anTestData, nNumNbrs, test_equal_training_data)


    if test_equal_training_data
        nNumNbrs_actual = nNumNbrs + 1;
    else 
        nNumNbrs_actual = nNumNbrs;
    end

    
    [nTestDim nTestSamples] = size(anTestData);
    [nTrainingDim nTrainingSample] = size(anTrainingData);

    anTestData_nbr_index = zeros(nNumNbrs, nTestSamples);
    anTestData_nbr_dist = zeros(nNumNbrs, nTestSamples);


    for i=1:nTestSamples


        [nbr_index, nbr_dist] = knnsearch(anTrainingData', anTestData(:,i)','k', nNumNbrs_actual, 'distance', 'minkowski');

        % if test data equals training data then
        % remove the first neighbor because it should be itself
        if test_equal_training_data
            nbr_index = nbr_index(2:length(nbr_index));
            nbr_dist = nbr_dist(2:length(nbr_dist));
        end


        anTestData_nbr_index(:,i) = nbr_index;
        anTestData_nbr_dist(:,i) = nbr_dist;
        
    end 

end




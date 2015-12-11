


function [anTestData_classify] = BAM_assign_knn_dist_label(anTestData_label_prob, ...
        anTestData_nbr_index, anTestData_nbr_dist, anTestData_nbr_dist_wt, DefaultLabelIndex)
    
% anTestData_label_prob
% anTestData_nbr_index
% anTestData_nbr_dist
% anTestData_nbr_dist_wt
% DefaultLabelIndex


    [nNumLabel nTestSamples] = size(anTestData_label_prob);
    

    anTestData_classify = zeros(1,nTestSamples);

    for i=1:nTestSamples
        
        [nMaxVal nMaxIndex] = max(anTestData_label_prob(:,i));
        
        if nMaxVal == 0
            nMaxIndex = DefaultLabelIndex;
        end
        
        anTestData_classify(i) = nMaxIndex;
               

    end


end




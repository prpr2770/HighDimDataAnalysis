

%% BAM_make_confusion_matrix
function [anConfusionMatrix, anPercentCorrect, nAllPercentCorrect] = ...
    BAM_make_confusion_matrix(anLabel, anTrueLabel, nNumLabel)

    %% assumes the labels are only from 1 to nNumLabel


    %% make confusion matrix
    anConfusionMatrix = zeros(nNumLabel, nNumLabel);
%    anPercentCorrect = zeros(1, nNumLabel);
%    nAllPercentCorrect = 0;

    nTotalSamples = length(anLabel);
    anLabelSamples = zeros(1, nNumLabel);
    
    % the number of samples given per label
    for i=1:nNumLabel
       anLabelSamples(i) = sum(anTrueLabel == i);
    end
    
    % making the confusion matrix
    for i=1:nTotalSamples

        nTrueIndex = anTrueLabel(i);
        nAssignedIndex = anLabel(i);
        
        anConfusionMatrix(nAssignedIndex, nTrueIndex) = ...
         anConfusionMatrix(nAssignedIndex, nTrueIndex) + 1;

    end

    anPercentCorrect = diag(anConfusionMatrix)' ./ anLabelSamples;
    nAllPercentCorrect = sum(diag(anConfusionMatrix)) / nTotalSamples;   

end


    
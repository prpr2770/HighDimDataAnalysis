function mutualSimilarity = getSimilarityMatrix(mutualDistance, similarityType, controlParameter);
% Define Similarity based on any function of mutualDistance matrix

% similarityType = 'InverseDistance'; % 'GaussianRBF' 'InverseDistance'

switch similarityType
    case 'InverseDistance'
        % Similarity = 1/(1 + Distance.^p)
        p = controlParameter;
        I = ones(size(mutualDistance));
        temp = I + mutualDistance.^p;
        mutualSimilarity = temp.^(-1);

    case 'GaussianRBF'
        sigmaSqr = controlParameter; %4*pi^2
        mutualSimilarity = exp(-mutualDistance.^2 ./ sigmaSqr); 
        
end

end

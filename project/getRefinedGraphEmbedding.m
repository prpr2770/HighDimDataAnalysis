function embeddedPoints = getRefinedGraphEmbedding(Dataset, dimReducedSpace)
%{
Derive Refined Graph Embedding: Using Gaussian RBF for Similarity
Computation.

%}

[totalPoints, dimDataPoints] = size(Dataset);
sigmaSqr = 4*pi^2;  % Gaussian RBF


% Determine Distances between the points and the nearest neighbors
mutualDistance = zeros(totalPoints,totalPoints);
mutualSimilarity = zeros(totalPoints,totalPoints);

for i=1:totalPoints
    for j=i:totalPoints
        dist_ij = sqrt(sum((Dataset(i,:) - Dataset(j,:)).^2));
        weight_ij = exp(-dist_ij^2/sigmaSqr);
        if i~=j
            % generate the distaneMatrix
            mutualDistance(i,j) = dist_ij;
            mutualDistance(j,i) = dist_ij;
            
            %generate the WeightMatrix
            mutualSimilarity(i,j) = weight_ij;
            mutualSimilarity(j,i) = weight_ij;
            
        else
            % diagonal entries
            mutualDistance(i,i) = dist_ij;
            mutualSimilarity(i,i) = weight_ij;
            
        end
        

    end
end

% ===========================================================
% Graph Theoretic Operations: Adjacency, Weight, Laplacian

Adj_M = mutualSimilarity;  % DETERMINE A THRESHOLDING FOR THIS MATRIX ???
deg_V = sum(Adj_M,2); %rowsum
Deg_M = diag(deg_V);

Laplacian_M = Deg_M - Adj_M;
normLaplacian_M = Deg_M^(-0.5)*(eye(size(Adj_M)) - Adj_M)*Deg_M^(-0.5);

probTrans_M = inv(Deg_M)*Adj_M;

% ===========================================================
% Graph Theoretic Dimension Reduction:

% dimReducedSpace = 3; % Dimension of Reduced Space

% Determine First 4 eigenvectors and eigenvalues of  normLaplacian
[eigVecs,Lambda_M] = eigs(normLaplacian_M,(dimReducedSpace+1));
lambda = diag(Lambda_M);

% Determine Stationary Distribution of Markov Chain: Is this a ROW VEC?
PI_vec = stationaryDist_MarkovChain(probTrans_M);
% extract required components of pi_Vec
pi_red_vec = PI_vec(1: dimReducedSpace + 1);



% Compute the fraction (1/sqrt(lambda.*pi_vec)
factor = lambda.*pi_red_vec;           % Column Vector
scaleFactor = 1./sqrt(factor(2:end)'); % Discard First Value; Row Vector
scaleFactor_M = repmat(scaleFactor,totalPoints,1);

% Compute: REFINED EMBEDDING
embeddedPoints = scaleFactor_M.*eigVecs(:,2:end);

end
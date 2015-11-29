function embeddedPoints = refinedGraphEmbedding_distanceMatrix(DistanceMatrix, dimReducedSpace)

totalPoints = length(DistanceMatrix);
mutualDistance = DistanceMatrix;

% mutualSimilarity = getSimilarityMatrix(mutualDistance, similarityType, controlParameter);
similarityType = 'GaussianRBF'; % ('GaussianRBF',{10, 100, 1000,...}) ('InverseDistance' , {1,2,3,...})
controlParam = 4*pi^2;
mutualSimilarity = getSimilarityMatrix(mutualDistance,similarityType, controlParam);

% ===========================================================
% Graph Theoretic Operations: Adjacency, Weight, Laplacian

Adj_M = mutualSimilarity;  % DETERMINE A THRESHOLDING FOR THIS MATRIX ???
deg_V = sum(Adj_M,2); %rowsum
Deg_M = diag(deg_V);

warning('computing Laplacian')
Laplacian_M = Deg_M - Adj_M;
normLaplacian_M = Deg_M^(-0.5)*(eye(size(Adj_M)) - Adj_M)*Deg_M^(-0.5);

probTrans_M = inv(Deg_M)*Adj_M;

% ===========================================================
% Graph Theoretic Dimension Reduction:

% dimReducedSpace = 3; % Dimension of Reduced Space

% Determine First 4 eigenvectors and eigenvalues of  normLaplacian
warning('computing eigenvectors')
[eigVecs,Lambda_M] = eigs(normLaplacian_M,(dimReducedSpace+1));
lambda = diag(Lambda_M); % eigenvalues

% Determine Stationary Distribution of Markov Chain: Is this a ROW VEC?
warning('computing MarkovChain stationary')
PI_vec = stationaryDist_MC(probTrans_M)
% extract required components of pi_Vec
pi_red_vec = PI_vec(1: dimReducedSpace + 1);


warning('computing Final Projection')
% Compute the fraction (1/sqrt(lambda.*pi_vec)
factor = abs(lambda).*pi_red_vec;            % Column Vector: ensure ABS value of Eigenvalues!
scaleFactor = 1./sqrt(factor(2:end)')   % Discard First Value; Row Vector
scaleFactor_M = repmat(scaleFactor,totalPoints,1);

% Compute: REFINED EMBEDDING
embeddedPoints = scaleFactor_M.*eigVecs(:,2:end);

end
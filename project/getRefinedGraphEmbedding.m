function embeddedPoints = getRefinedGraphEmbedding(Dataset, dimReducedSpace)
%{
Derive Refined Graph Embedding: Using Gaussian RBF for Similarity
Computation.

%}

[totalPoints, dimDataPoints] = size(Dataset);
sigmaSqr = 4*pi^2;  % Gaussian RBF
epsilon = 0;

% Determine Distances between the points and the nearest neighbors
mutualDistance = zeros(totalPoints,totalPoints);
mutualSimilarity = zeros(totalPoints,totalPoints);

tic
for i=1:totalPoints
    i
    ref_vec = Dataset(i,:);
    ref_vec_repMat = repmat(ref_vec,totalPoints,1);
    
    try
        dist_ijs = (sum((ref_vec_repMat - Dataset).^2, 2) ).^(0.5); %row_sum
        weight_ijs = exp(-dist_ijs.^2 / sigmaSqr);
    catch 
        warning('|ref_vec_repMat| != |Dataset|')
    end

    % generate distanceMatrix
    mutualDistance(i,:) = dist_ijs';
    mutualDistance(:,i) = dist_ijs;
    
    % generate weightMatrix
    mutualSimilarity(i,:) = weight_ijs';
    mutualSimilarity(:,i) = weight_ijs;
    
    % -----------------------------------------------------------

    
end
toc
% ===========================================================
% Graph Theoretic Operations: Adjacency, Weight, Laplacian


Adj_M = ones(size(mutualSimilarity)).*(mutualSimilarity > epsilon);  % DETERMINE A THRESHOLDING FOR THIS MATRIX ???
deg_V = sum(Adj_M ,2); %rowsum
Deg_M = diag(deg_V);

warning('Laplacian')
Laplacian_M = Deg_M - Adj_M;
normLaplacian_M = Deg_M^(-0.5)*(eye(size(Adj_M)) - Adj_M)*Deg_M^(-0.5);

probTrans_M = inv(Deg_M)*Adj_M;

% ===========================================================
% Graph Theoretic Dimension Reduction:

% dimReducedSpace = 3; % Dimension of Reduced Space

% Determine First 4 eigenvectors and eigenvalues of  normLaplacian
[eigVecs,Lambda_M] = eigs(normLaplacian_M,(dimReducedSpace+1));
lambda = diag(Lambda_M);

warning('Stationary Markov Chain')
tic
% Determine Stationary Distribution of Markov Chain: Is this a ROW VEC?
% PI_vec = stationaryDist_MarkovChain(probTrans_M);
PI_vec = stationaryDist_MC(probTrans_M);
toc
% extract required components of pi_Vec
pi_red_vec = PI_vec(1: dimReducedSpace + 1);



% Compute the fraction (1/sqrt(lambda.*pi_vec)
factor = lambda.*pi_red_vec;           % Column Vector
scaleFactor = 1./sqrt(factor(2:end)'); % Discard First Value; Row Vector
scaleFactor_M = repmat(scaleFactor,totalPoints,1);

warning('Refined Embedding')
% Compute: REFINED EMBEDDING
embeddedPoints = scaleFactor_M.*eigVecs(:,2:end);

end
% Test Cases: 
% ---------------------------------
% 
% A)Two intersecting Spherical Shells: Using Graph embedding to separate them out!
% B)Two intersecting Spheres: Using Graph embedding to separate them out!

clear all; close all; clc;

testScenario = 'A'; %'B' 

N = 10000; % Number of points in 1 cluster
r = 1;
n = 1000; %dimensions of space
alpha = 1.7;    % fraction of diameter of sphere to shift the second sphere

% ===========================================================
% Generate the Dataset

switch testScenario
    case 'A'
        % Case A: Intersecting Spherical Shells
        %--------------------------------------------------
        validSamplesOnBall_1 = samplePointsInBall_Uniform(r,n,N); % Shell 1
        validSamplesOnBall_2 = samplePointsInBall_Uniform(r,n,N); % Shell 2

    case 'B'
        % Case B: Intersecting Spheres
        %--------------------------------------------------
        validSamplesInBall_1 = samplePointsInBall_UniformSolid(r,n,N);
        validSamplesInBall_2 = samplePointsInBall_UniformSolid(r,n,N);
        
end

% Displace second shell by a distace 1.8*r
dirVector = rand(1,n);
dirVector = alpha*(2*r)* dirVector / norm(dirVector,2);
dirVector_Matrix = repmat(dirVector,N,1);
validSamplesOnBall_2 = validSamplesOnBall_2 +   dirVector_Matrix;
    
% Create a complete dataset from the two
Dataset = [validSamplesOnBall_1; validSamplesOnBall_2];
% ===========================================================
% Refined Graph Embedding: Using Gaussian RBF
dimReducedSpace = 3;
graphEmbeddedPoints = getRefinedGraphEmbedding(Dataset, dimReducedSpace)

% ===========================================================
% Derive FJLT based Dimension Reduction


% ===========================================================
% Generate Plots:

g_X = graphEmbeddedPoints(:,1);
g_Y = graphEmbeddedPoints(:,2);
g_Z = graphEmbeddedPoints(:,3);

scatter3(g_X, g_Y, g_Z)



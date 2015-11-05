% Author: Prasanth Prahladan
% HW1: Assignment 3.3
% Algorithm 2:Generation of Samples with Uniform Measure on S_(n-1)(1) with
% using GAUSSIAN MEASURE

function validSamplesOnBall = samplePointsInBall_Uniform(r,n,N)
% n:                    number of dimensions
% N:                    number of samples
% r:                    radius of Ball
% validSamplesInBall:   list of valid samples obtained via Rejection Sampling                  

g = randn(N,n);                     % gaussian distributed points - row vectors
sampleNorms = sqrt(sum(g.^2,2));    % rowSum = 2

% repeat norms for each dimension, to obtain Nxn matrix
normsMat = repmat(sampleNorms,1,n);

% compute the Points on Surface of Ball of Radius R
validSamplesOnBall = r.*g./normsMat;



end
% Author: Prasanth Prahladan
% HW1: Assignment 4
% Algorithm 3: Uniform sampling in Ball 

function validSamplesInBall = samplePointsInBall_UniformSolid(r,n,N)
% n:                    number of dimensions
% N:                    number of samples
% r:                    radius of Ball
% validSamplesInBall:   list of valid samples obtained via Rejection Sampling                  

g = randn(N,n);                     % gaussian distributed points
u = rand(N,1);
u = u.^(1/n);

sampleNorms = sqrt(sum(g.^2,2));    % rowSum = 2

% repeat norms for each dimension, to obtain Nxn matrix; 
% This facilitates division to obtain unit-vector.
normsMat = repmat(sampleNorms,1,n);

% repeat the uniform distribution samples
uniMat = repmat(u,1,n);

% compute the Points on Surface of Ball of Radius R
validSamplesInBall = r.*uniMat.*g./normsMat;

end
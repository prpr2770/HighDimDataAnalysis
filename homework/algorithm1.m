% Author: Prasanth Prahladan
% HW 1: Assignment 2
% 1. Generate 10000 samples in unit-cube of n-dim space: n = {1, ..., 400}
% 2. Plot number of points rejected, i.e not lying inside Unit-Ball, as
% n-varies.

N = 10000;
r = 1;
n = 20;

iterMax = 1000;
pointsRejected = zeros(iterMax,n);

% Q1: 
for iter=1:iterMax
iter
    for dim=1:n    
        [numPointsRejected, validSamplesInBall] = samplePointsInBall(r,dim,N);
        pointsRejected(iter,dim) = numPointsRejected;
    end

end
n_vary = 1:1:n;

avgPointsRejected = sum(pointsRejected,1)/iterMax;

% Q2:
plot(n_vary,avgPointsRejected)
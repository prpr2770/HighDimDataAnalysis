% Author: Prasanth Prahladan
% HW 1: Assignment 2
% Algorithm 1:Generation of Samples with Uniform Measure in B_n(1) with
% REJECTION METHOD: UNIFORM MEASURE

function [numPointsRejected, validSamplesInBall] = samplePointsInBall(r,n,N)
% n:                    number of dimensions
% N:                    number of samples
% r:                    radius of Ball
% numPointsRejected:    number of points rejected
% validSamplesInBall:   list of valid samples obtained via Rejection Sampling                  

numPointsRejected = 0;
countPointsInBall = 0;

for i=1:N
   u = rand(1,n);
   x = r*(1-2*u);
   if (norm(x,2) > r)
       numPointsRejected = numPointsRejected + 1;
       validSamplesInBall(1,:) = zeros(1,n);
   else
       countPointsInBall = countPointsInBall + 1;
       validSamplesInBall(countPointsInBall,:) = x;
   end
    
end
% Author: Prasanth Prahladan
% HW 1: Assignment 3.1
% Generation of Samples with GAUSSIAN MEASURE in R_N.

function validSamples = sampleGaussianPoints(n, mean, stdev, N)
% n:                    number of dimensions: SQUARED INTEGER < 400
% N:                    number of samples
% validSamples:         list of valid samples 

    if sqrt(n) ~= floor(sqrt(n))
        fprintf('Error! n should be Squared Integer!\n');
        return
    else
        x = mean + stdev.*randn(N,n);
        validSamples = x;
    end

end


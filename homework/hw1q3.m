% Author: Prasanth Prahladan
% HW 1: Assignment 3.1
% 3. Generate 10000 samples x in R_n from Gaussian Distribution, 0 Mean,
% unit Var, with n = squared_integer <400
% 4. Plot histogram of norm(x) for n = {4,25,100, 225, 400}
% 5. Compute Mean, Variance of each distribution and plot it as a function
% of n. What do you notice? Make a conjecture on the concentration of
% norm(x), when x is gaussian distributed in R_n


N = 10000;
n_vals = [4, 25, 100, 225, 400];
%n_vals = 5;

% Generate Gaussian Distribution
mean = 0;
stdev = 1;

% Archiving variables
listSampleNorms = zeros(N,length(n_vals));
meanSamples = zeros(1,length(n_vals));
varSamples  = zeros(1,length(n_vals));


for i=1:length(n_vals)    
    n = n_vals(i);
    validSamples = sampleGaussianPoints(n, mean, stdev, N);

    % compute 2-norms of the samples
    sampleNorms = sqrt(sum(validSamples.^2,2));    % rowSum = 2
    listSampleNorms(:,i) = sampleNorms;
   
    % compute Mean
    meanSamples(i) = sum(sampleNorms,1)/N;                % colSum = 1
    
    % compute Variance: E[X^2] - E[X]^2
%     varSamples(i) = sum(sampleNorms.^2,1)/N - meanSamples(i)^2;
    varSamples(i) = var(sampleNorms);
end

fig1 = figure(1)
subplot(2,1,1)
plot(n_vals, meanSamples,'b--o');
ylabel('Mean')
xlabel('Dimension')
subplot(2,1,2)
plot(n_vals, varSamples,'r-*');
ylabel('Variance')
xlabel('Dimension')

fig2 = figure(2)
[H,X] = hist(listSampleNorms,100);
plot(X, H)
% Figure properties
str = num2str(n_vals);
str = strsplit(str,' ');
legend(str)
xlabel('Norm Magnitude')
ylabel('Frequency')


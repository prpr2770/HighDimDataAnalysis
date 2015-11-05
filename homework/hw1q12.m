% Author: Prasanth Prahladan
% HW 1: Assignment 3.3
% 12. Generate 10000 points on Sphere(S_(n-1)(1)) for n={1,..,400}. Project
% n points on axis(x1). 
% Plot histogram of the projections for n = {4,25,100, 225, 400}.
% 13. Compute mean and variance of the distribution of the projections, and
% plot these values as a function of n.


N = 10000; %total number of points
r = 1;
n = 400; %dimensions of space



% --------------------------------------------------------------
n_vals = [4, 25, 100, 225, 400];

listSampleProjs = zeros(N,length(n_vals));
meanSamples = zeros(1,length(n_vals));
varSamples = zeros(1,length(n_vals));


for i = 1:length(n_vals)
    n = n_vals(i);
    validSamplesOnBall = samplePointsInBall_Uniform(r,n,N);
    
    % project samples on axis x1
    samples_proj_on_X1 = validSamplesOnBall(:,1);
    
    % Archive data + stats
    listSampleProjs(:,i) = samples_proj_on_X1;
    meanSamples(i) = sum(samples_proj_on_X1)/N;
    varSamples(i) = sum(samples_proj_on_X1.^2)/N - meanSamples(i)^2;

end
% --------------------------------------------------------------
% Plots:

fig1 = figure(1)
subplot(2,1,1)
plot(n_vals, meanSamples,'b--o');
ylabel('Mean')
xlabel('Dimension')
subplot(2,1,2)
plot(n_vals, varSamples,'r-*');
ylabel('Variance')  
xlabel('Dimension')
title('Mean(Average) and Variance of the Samples')



fig2 = figure(2)
numBins = 1000;
[H,X] = hist(listSampleProjs,numBins);
plot(X, H)
% Figure properties
str = num2str(n_vals);
str = strsplit(str,' ');
legend(str)
xlabel('Projected Length ')
ylabel('Frequency')
title('Histogram of Samples for different dimensions')


fig3 = figure(3)
plot(X,'--o')
title('Histogram Bin centers')


fig5 = figure(5)
posEpsVals = zeros(1,length(n_vals));
negEpsVals = zeros(1,length(n_vals));
for ind = 1:length(n_vals)
    total = zeros(1,numBins/2);
    xend = zeros(1,numBins/2);
    for j=1:numBins/2;
       x = X(j:numBins-j+1);
       h = H(j:numBins-j+1,ind);
       total(j) = sum(h);
       xend(j) = x(end);
       if total(j) >= 0.99*N
           posEpsVals(ind) = x(end);
           negEpsVals(ind) = x(1);
       end
    end
    plot(xend,total)
    hold on
end
plot(xend,0.99*N*ones(size(total)),'r')
hold off
xlabel('Eps Bandwidth')
ylabel('Sample Counts')
title('Total datasamples contained within bandwidth')


fig4 = figure(4)
plot(n_vals,posEpsVals,'g--*')
hold on
plot(n_vals,negEpsVals,'r--*')
hold off
title('Upper and Lower Eps Bounds to contain 0.99N samples')
    

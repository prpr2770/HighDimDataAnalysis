% Author: Prasanth Prahladan
% HW 1: Assignment 4
% 16. Generate N = 10000 points inside Ball_n(sqrt(n)) for n={1,..,400}
% 17. Project points onto x1 coordinate. Plot histogram for n={4,25,100,225,400}
% 18. Compute Mean and Variance of Distribution of the projections, and plot these wrt n
% 19. 


N = 10000;
numBins = 100;

n_vals = [4, 25, 100, 225, 400];

listSampleProjs = zeros(N,length(n_vals));
meanProjSamples = zeros(1,length(n_vals));
varProjSamples = zeros(1,length(n_vals));

listSampleNorms = zeros(N,length(n_vals));


for i = 1:length(n_vals)
    n = n_vals(i);
    r = sqrt(n);
    validSamplesInBall = samplePointsInBall_UniformSolid(r,n,N);

    % ----------------------------------------------------
    % project samples on axis x1
    proj_samples_on_X1 = validSamplesInBall(:,1);
    % Archive data + stats
    listSampleProjs(:,i) = proj_samples_on_X1;
    % Computing mean and variance
    meanProjSamples(i) = mean(proj_samples_on_X1);
    varProjSamples(i) = var(proj_samples_on_X1);

    % ----------------------------------------------------
    % Similar computations as above for Sample-Norms!
    norm_samples_on_X1 = sqrt(sum(validSamplesInBall.^2 ,2));
    listSampleNorms(:,i) = norm_samples_on_X1;
    meanNormSamples(i) = mean(norm_samples_on_X1);
    varNormSamples(i) = var(norm_samples_on_X1);
end
% --------------------------------------------------------------
% --------------------------------------------------------------
% Plots: for Projections of Samples onto X1

fig1 = figure(1)
subplot(2,1,1)
plot(n_vals, meanProjSamples,'b--o');
ylabel('Mean')
xlabel('Dimension')
subplot(2,1,2)
plot(n_vals, varProjSamples,'r-*');
ylabel('Variance')  
xlabel('Dimension')
title('Mean(Average) and Variance of the Sample-Projection')


fig2 = figure(2)
[H,X] = hist(listSampleProjs,numBins);
plot(X, H)
% Figure properties
str = num2str(n_vals);
str = strsplit(str,' ');
legend(str)
xlabel('Projected Length ')
ylabel('Frequency')
title('Histogram of Samples for different dimensions')

% Q:19:
fig3 = figure(3)
geqLowerLimit = ( listSampleProjs >= -0.5 );
leqUpperLimit = ( listSampleProjs <= 0.5 );
valuesInRange = geqLowerLimit.*leqUpperLimit;
relVolumes = sum(valuesInRange,1)/N;

plot(n_vals, relVolumes);
ylabel('Relative Volume')
xlabel('#Dimensions')
%title('Plot of relative Volumes of B_n( $$\sqrt{n}$$ )','interpreter','latex')
title('Plot of relative Volumes of $$B^n(\surd{n})$$','interpreter','latex' )


% --------------------------------------------------------------
% --------------------------------------------------------------
% Plots for Norms of Samples


fig11 = figure(11)
subplot(2,1,1)
plot(n_vals, meanNormSamples,'b--o');
ylabel('Mean')
xlabel('Dimension')
subplot(2,1,2)
plot(n_vals, varNormSamples,'r-*');
ylabel('Variance')  
xlabel('Dimension')
title('Mean(Average) and Variance of the Sample-Norm')


fig12 = figure(12)
numBins = 1000;
[H,X] = hist(listSampleNorms,numBins);
plot(X, H)
% Figure properties
str = num2str(n_vals);
str = strsplit(str,' ');
legend(str)
xlabel('Norm ')
ylabel('Frequency')
title('Histogram of Sample-Norms for different dimensions')

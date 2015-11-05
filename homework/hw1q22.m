% HW 1: 
% Prasanth Prahldan
% Q.22/23 Create Random Wigner Matrices
% a) Symmetric Bernoulli Ensemble with p=0.5
% b) Gaussian Orthogonal Ensembles(GOE)

clear all;close all;clc;

TotalBins = 20;
TotalSamples = 100;
% n_vals = [10, 50, 100, 500, 1000];
n_vals = [10, 50];


WG_minL = zeros(1,length(n_vals));
WG_maxL = zeros(1,length(n_vals));

BE_minL = zeros(1,length(n_vals));
BE_maxL = zeros(1,length(n_vals));


for iter = 1: length(n_vals)
   n = n_vals(iter)
   

   % ---------------------------------------
   % For realizations
   WG_lambdaMin = zeros(1,TotalSamples);
   WG_lambdaMax = zeros(1,TotalSamples);
   
   BE_lambdaMin = zeros(1,TotalSamples);
   BE_lambdaMax = zeros(1,TotalSamples);
   
   BE_NrmEigs = zeros(n,TotalSamples);
   WG_NrmEigs = zeros(n,TotalSamples);
   
   for sample = 1:TotalSamples

       % obtain Matrices
       BE = generateSymmBernoulliEnsemble(n,0.5);
       WG = generateGaussianOrthoEnsemble(n);

       % compute the eigenvalues: 'BE' = Both-Ends
       BE_MinMaxL = eigs(BE,2,'BE');
       WG_MinMaxL = eigs(WG,2,'BE');

       % -----------------------------------
       % Archiving
       WG_lambdaMin(sample) = WG_MinMaxL(1);
       WG_lambdaMax(sample) = WG_MinMaxL(2);
       
       BE_lambdaMin(sample) = BE_MinMaxL(1);
       BE_lambdaMax(sample) = BE_MinMaxL(2);
       
      % -----------------------------------
      % Computing the ESD - Derive the eigenvalues 
      BE_NrmEigs(:,sample) = eig(BE)/sqrt(n);
      WG_NrmEigs(:,sample) = eig(WG)/sqrt(n);
       
   end
   
   figure(iter*11)
   hist(BE_NrmEigs,100)
   title('Histogram of normalized eigenvalues for n=%d',n)

   figure(iter*13)
   hist(WG_NrmEigs,100)
   title('Histogram of normalized eigenvalues for n=%d',n)
   
   % Find min and Max values for each n
   WG_minL(iter) = min(WG_lambdaMin);
   WG_maxL(iter) = max(WG_lambdaMax);
   
   BE_minL(iter) = min(BE_lambdaMin);
   BE_maxL(iter) = max(BE_lambdaMax);
   
end

% -------------------------------------------------------------
% Determine Polynomial fit values for the eigenvalues
polyorder = 4; 

% Wigner Gaussian
% Compute polynomial coefficients
pf_WG_minL = polyfit(n_vals,WG_minL,polyorder);
pf_WG_maxL = polyfit(n_vals,WG_maxL,polyorder);
% Compute the polynomial approximators
xf_WG_minL = polyval(pf_WG_minL , n_vals);
xf_WG_maxL = polyval(pf_WG_maxL, n_vals);

% Wigner Bernoulli
% Compute polynomial coefficients
pf_BE_minL = polyfit(n_vals,BE_minL,polyorder);
pf_BE_maxL = polyfit(n_vals,BE_maxL,polyorder);
% Compute the polynomial approximators
xf_BE_minL = polyval(pf_BE_minL , n_vals);
xf_BE_maxL = polyval(pf_BE_maxL, n_vals);

% -------------------------------------------------------------
fig1 = figure(1)
scatter(n_vals,WG_minL,'bo')
hold on
plot(n_vals, xf_WG_minL,'b--')
hold on
scatter(n_vals,WG_maxL,'ko')
hold on
plot(n_vals, xf_WG_maxL,'k--')
hold off
title('Min-Max EVs for Wigner Gaussian')
xlabel('Spatial Dimension n')

fig2 = figure(2)
scatter(n_vals,BE_minL,'bo')
hold on
plot(n_vals, xf_BE_minL,'b--')
hold on
scatter(n_vals,BE_maxL,'ko')
hold on
plot(n_vals, xf_BE_maxL,'k--')
hold off
title('Min-Max EVs for Bernoulli')
xlabel('Spatial Dimension n')


% -------------------------------------------------------------
% Q24:Line of Best-Fit: Logarithm of n - vs - Eigenvalues Plot
% -------------------------------------------------------------
% Determine Polynomial fit values for the eigenvalues
polyorder = 1; 

log_n_vals = log(n_vals);
% Wigner Gaussian
% Compute polynomial coefficients
pf_WG_minL = polyfit(log_n_vals,WG_minL,polyorder)
pf_WG_maxL = polyfit(log_n_vals,WG_maxL,polyorder)
% Compute the polynomial approximators
xf_WG_minL = polyval(pf_WG_minL , log_n_vals);
xf_WG_maxL = polyval(pf_WG_maxL, log_n_vals);

% Wigner Bernoulli
% Compute polynomial coefficients
pf_BE_minL = polyfit(log_n_vals,BE_minL,polyorder)
pf_BE_maxL = polyfit(log_n_vals,BE_maxL,polyorder)
% Compute the polynomial approximators
xf_BE_minL = polyval(pf_BE_minL , log_n_vals);
xf_BE_maxL = polyval(pf_BE_maxL, log_n_vals);

% -------------------------------------------------------------
fig11 = figure(11)
scatter(log_n_vals,WG_minL,'bo')
hold on
plot(log_n_vals, xf_WG_minL,'b--')
hold on
scatter(log_n_vals,WG_maxL,'ko')
hold on
plot(log_n_vals, xf_WG_maxL,'k--')
hold off
title('Min-Max EVs for Wigner Gaussian')
xlabel(' log(n) ')



fig12 = figure(12)
scatter(log_n_vals,BE_minL,'bo')
hold on
plot(log_n_vals, xf_BE_minL,'b--')
hold on
scatter(log_n_vals,BE_maxL,'ko')
hold on
plot(log_n_vals, xf_BE_maxL,'k--')
hold off
title('Min-Max EVs for Bernoulli')
xlabel(' log(n)')



fig13 = figure(13)
scatter(log_n_vals,WG_minL,'bo')
hold on
plot(log_n_vals, xf_WG_minL,'b--')
hold on
scatter(log_n_vals,WG_maxL,'bo')
hold on
plot(log_n_vals, xf_WG_maxL,'b--')
hold on
scatter(log_n_vals,BE_minL,'k+')
hold on
plot(log_n_vals, xf_BE_minL,'k-.')
hold on
scatter(log_n_vals,BE_maxL,'k+')
hold on
plot(log_n_vals, xf_BE_maxL,'k-.')
hold off
title('Min-Max EVs for Wigner Bernoulli and Gaussian')
xlabel(' log(n) ')


% -------------------------------------------------------------
% -------------------------------------------------------------






fig3 = figure(3)
x = 0:0.1:3;
mask = abs(x)<= 2;
esd = 1/(2*pi)* (4 - x.^2).^(0.5).*mask;
plot(x,esd,'r')
title('Semi-Circular Distribution')

%{
% Q25: Compute Empirical Spectral Distribution

fig4 = figure(4)
subplot(2,1,1)
[P,X] = hist(nrmLambda_BE_n,TotalBins);
P = P/length(nrmLambda_BE_n);
plot(X,P)
hold on
plot(x,esd,'r')
title('ESD: Bernoulli Matrix')

subplot(2,1,2)
[P,X] = hist(nrmLambda_WG_n,TotalBins);
P = P/length(nrmLambda_WG_n);
plot(X,P)
hold on
plot(x,esd,'r')
title('ESD: Gaussian Matrix')

%}
% HW 1: 
% Prasanth Prahldan

% Q.26. Display normalized eigenvalues using ScatterPlot
% Generate for Random NON-SYMMETRIC matrices!

clear all;close all;clc;

TotalBins = 20;
TotalSamples = 100;
%n_vals = [10, 50, 100, 500, 1000];
n_vals = [10, 25, 50];

   % list of all eigenvalues for all n
   nrmLambda_BR_n = cell(1,length(n_vals)); % Cell-Array
   nrmLambda_GR_n = cell(1,length(n_vals)); % Cell-Array


for iter = 1: length(n_vals)
   n = n_vals(iter)


   
   % list of all eigenvalues for all realizations
   nrmLambda_BR = [];
   nrmLambda_GR = [];
   
   for sample = 1:TotalSamples

       % obtain Matrices
       BR = generateBernoulliRandomMatrix(n,0.5);
       GR = generateGaussianRandomMatrix(n);

       % normalized eigenvalues: Column Vectors
       smpl_nrmLambda_BR = 1/sqrt(n)*eig(BR);
       smpl_nrmLambda_GR = 1/sqrt(n)*eig(GR);
       
       % compile list of all normalized eigenvalues
       nrmLambda_BR = [nrmLambda_BR; smpl_nrmLambda_BR];
       nrmLambda_GR = [nrmLambda_GR; smpl_nrmLambda_GR];
       
   end
   
   % Compile the normalized evs for each n: They are of different sizes
   nrmLambda_BR_n{iter} = nrmLambda_BR; 
   nrmLambda_GR_n{iter} = nrmLambda_GR;
end


colors = ['g','r','b','k','y','cyan','magenta'];
%Eigen-values are complex. So how would this plot out?
%For each n, create different scatter-plots - on Real-Vs-Imag components.

fig1 = figure(1)
for i=1:length(nrmLambda_BR_n)
    real_nrnLambda_BR_n = real(nrmLambda_BR_n{i});
    imag_nrnLambda_BR_n = imag(nrmLambda_BR_n{i});
    scatter(real_nrnLambda_BR_n,imag_nrnLambda_BR_n,colors(i),'+')
    hold on
end
hold off
title('Plot of eigenvalues for Bernoulli Random Matrix')

fig1a = figure(11)
subplot(2,1,1)
[Nr,Xr] =hist(real_nrnLambda_BR_n,100)
plot(Xr, Nr)
xlabel('Real Component')
subplot(2,1,2)
[Ni,Xi]= hist(imag_nrnLambda_BR_n,100)
plot(Xi, Ni)
xlabel('Imaginary Component')
title('Histogram of real and imag components of evs for Bernoulli Random Matrices')
fig2 = figure(2)

for i=1:length(nrmLambda_GR_n)
    real_nrnLambda_GR_n = real(nrmLambda_GR_n{i});
    imag_nrnLambda_GR_n = imag(nrmLambda_GR_n{i});
    scatter(real_nrnLambda_GR_n,imag_nrnLambda_GR_n,colors(i),'+')
    hold on
end
hold off
title('Plot of eigenvalues for Gaussian Random Matrix')

fig2a = figure(21)
subplot(2,1,1)
[Mr, Yr] = hist(real_nrnLambda_GR_n,100)
plot(Yr,Mr)
xlabel('Real Component')
subplot(2,1,2)
[Mi, Yi] = hist(imag_nrnLambda_GR_n,100)
plot(Yi,Mi)
xlabel('Imaginary Component')
title('Histogram of real and imag components of evs for Gaussian Random Matrices')
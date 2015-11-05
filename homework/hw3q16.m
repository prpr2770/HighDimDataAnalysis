%{
HW#3 Q16. 
For n “ 300, and for each combination of ? “ 5, 6, . . . , 50 and ? “ 1, 2, . . . , 50 generate 20
random realizations of the planted partition model. Use the algorithm Partition to detect the
communities, and compute the overlap score overlap.
You will represent the results of this experiment by constructing a matrix score(?, ?) that will
contain the mean (computed over the 20 experiments) of the overlap score.
Display the matrix in grayscale using the MATLAB function imagesc, and overlay the implicit
curve dened by (20) that denes the condition on ? on ? for the detectability of the communities
(see Fig. 6-left for an example of the gure).


Q17. 
For n “ 300, and for each combination of a “ 5, 6, . . . , 70 and b “ 1, 2, . . . , 50 generate 20
random realizations of the planted partition model. Use the algorithm Partition to detect the
communities, and compute the overlap score overlap.
You will represent the results of this experiment by constructing a matrix score(?, ?) that will
contain the mean (computed over the 20 experiments) of the overlap score.

%}
close all; clear all; clc
tic
n = 300;
numTrials = 20;
alphaMax = 70;
betaMax = 50;
alphaMin = 5;
betaMin = 1;

Alpha = alphaMin:1:alphaMax;
Beta = betaMin:1:betaMax;

overlapMatrix = zeros(length(Alpha), length(Beta));

for i=1:length(Alpha)
    for j= 1:length(Beta)
        
        alpha = Alpha(i);
        beta = Beta(j);
        
% %         % Dense Matrix
% %         p = alpha/n*log(n);
% %         q = beta/n*log(n);
         
        % Sparse Matrix
        p = alpha/n;
        q = beta/n;
        
        overlapScore = zeros(1,numTrials);
        % ---------------------------------------------------------
        % We need to run the algorithm only if q<p
        % We set the default values to zero!
        if (q<p)
         for iter = 1:numTrials
             [A, w] = getPartitionGraphModel(n,p,q);
             w_pred = runPartitionAlgo(A);
             overlapScore(iter) = getPartitionOverlap(w, w_pred);
         end 
        end
        % ---------------------------------------------------------
         % store the avg. Overlap score for given (alpha,beta)
         overlapMatrix(i,j) = sum(overlapScore)/ numTrials;
    end 
end

% % % % ----------------------------------------------------
% % % % Plot the following curve on the 
% % % % alpha - beta > sqrt(0.5*(alpha + beta))*2 / sqrt(log(n))
% % % % k = 2/sqrt(log n)
% % % % alpha > 0.5*( (2*beta + k^2/2)\pm k/2 sqrt(16*beta + k^2/2)   )
% % % 
% % % k = 2/sqrt(log(n));
% % % Alpha1 = 0.5*((2*Beta + k^2/2) + (k/2)*sqrt(16*Beta + k^2));
% % % Alpha2 = 0.5*((2*Beta + k^2/2) - (k/2)*sqrt(16*Beta + k^2));
% % % % ----------------------------------------------------

% ----------------------------------------------------
% How does this change for the sparse matrix?
% SPARSE MATRIX BOUNDARIES
Alpha1 = (Beta + 1) + sqrt(1+ 4*Beta);
Alpha2 = (Beta + 1) - sqrt(1+ 4*Beta);

% ----------------------------------------------------
close all

fig1 = figure(1)
% Plot the image as a GRAYSCALE
betaDim = [betaMin betaMax];
alphaDim = [alphaMin alphaMax];
img = imagesc(betaDim, alphaDim, overlapMatrix);
colormap(gray);
colorbar;
hold on
ax2 = plot(Beta,Alpha1,'r','LineWidth',2);
% colormap(ax2, parula )
hold on
ax3 = plot(Beta,Alpha2,'r','LineWidth',2);
% colormap(ax3, parula )
hold off
title('Sparse Matrix: Probability of successfully detecting Partitions')
% title('Dense Matrix: Probability of successfully detecting Partitions')


fig2 = figure(2)
% Plot the image as a GRAYSCALE
betaDim = [betaMin betaMax];
alphaDim = [alphaMin alphaMax];
img = imagesc(betaDim, alphaDim, overlapMatrix);
colormap(gray);
colorbar;
hold on
ax2 = plot(Beta,Alpha1,'r','LineWidth',3);
% colormap(ax2, parula )
hold on
ax3 = plot(Beta,Alpha2,'r','LineWidth',3);
% colormap(ax3, parula )
hold off
title('Sparse Matrix: Probability of successfully detecting Partitions')
% title('Dense Matrix: Probability of successfully detecting Partitions')
toc
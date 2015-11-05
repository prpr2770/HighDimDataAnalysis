clear all;close all;clc

N=1000;
n = 5;

% ==============================================
% Generate Data Points:
% Vectors on Unit-Sphere in 5-D space
data = zeros(1,5,N);
for i=1:N
    A = rand(1,n);
    data(:,:,i) = A/sum(A);
end
% ==============================================
% Number of pairs being compared
M = N*(N-1)/2;

rsqr = zeros(1,M);
count = 0;
for i=1:N
for j=i+1:N
    count = count +1;
    dist = norm(data(:,:,i)-data(:,:,j),'fro');
    rsqr(count) = dist;
end
end

% ==============================================
rMin = min(rsqr);
rMax = max(rsqr);

% Determine RMAX -> Used for Scaling.
RMAX = ceil(rMax) + 1;
rsqr_Scaled = rsqr/RMAX;

% r-Values for Plotting
rVals = 0.01:0.01:0.4;
C = zeros(size(rVals));

for i = 1:length(rVals)
    r = rVals(i);
    C(i) = sum(rsqr_Scaled < r*ones(size(rsqr_Scaled)));
end

C = C/M;

dims = log(C)./log(rVals);
    
plot(rVals,2.^dims)
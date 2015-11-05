function X = GenerateSphere(n,k,D)

%
% function X = GenerateBall(n,k,D)
% 
% Generates n uniformly distributed points in the k-dimensional unit ball in R^D
%
% IN:
%   n   : number of points
%   k   : dimension of the ball
%   D   : ambient dimension
%
% OUT:
%   X   : n by D matrix of points
%

k=k+1;  % a k dimensional sphere exists in k+1 dimensional space


% Draw points on the unit sphere first (uniform angular)
X       = [randn(n,k),zeros(n,D-k)];
Xnorms  = sqrt(sum(X.^2,2));

for i=1:n
    X(i,:)=X(i,:)/Xnorms(i);
end

end
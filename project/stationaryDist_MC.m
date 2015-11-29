function [y2]=stationaryDist_MC(Q)
% Q is the transition matrix of the Markov Chain. We use the method of state
% space reduction to compute the stationary distributions of Q.

% Source: http://www.mathworks.com/matlabcentral/answers/uploaded_files/4107/stationary.m

syms z; 
P=Q;
[ns ms]=size(P);
n=ns;
n8=n;


while n8>1 % total states > 1 
    n8;
Q = P(1:n8-1, 1:n8-1);
colVec = P(1:n8-1,n8);
rowVec = P(n8,1:n8-1);

colMat = repmat(colVec,1,n8-1);
rowMat = repmat(rowVec,n8-1,1);
scaleFactor = 1-P(n8,n8);
scaleMat = 1/scaleFactor*ones(size(rowMat));

% compute the update
Q = Q + colMat.*rowMat.*scaleMat;

% reassign values
P(1:n8-1,1:n8-1) = Q;

n8=n8-1;
end

%{
% ------------------------------------------------------------
% Why is this being done symbollically?
y=sym(ones(n,1));
n3=1;
while n3<n
    n3
    x=sym(ones(n3,1));
    x(1:n3) = y(1:n3);
            
    z=sum(x.*P(1:n3,n3+1))/(sym(1)-P(n3+1,n3+1));
    y(n3+1)=z/(sym(1)+z);
    
    y(1:n3) = (sym(1)-y(n3+1))* y(1:n3);
    n3=n3+1;
end
%}

% ------------------------------------------------------------
y=ones(n,1);
n3=1;
while n3<n
    n3;
    x = y(1:n3);
    
%     size(x)
%     size(P(1:n3,n3+1))      
    
    z=sum(x.*P(1:n3,n3+1))/(1-P(n3+1,n3+1));
    y(n3+1)=z/(1+z);
    
    y(1:n3) = (1-y(n3+1))* y(1:n3);
    n3=n3+1;
end


y2=y;
end
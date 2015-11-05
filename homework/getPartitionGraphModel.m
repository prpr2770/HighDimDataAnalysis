% ECEN 5322: HW#3 - Experiments

function [A, partitionIndicatorVec] = getPartitionGraphModel(n,p,q)
%{ 
Q 13.
n : total #nodes in G
p : probability of link between two vertices inside Cluster
q : probability of link edges between two vertices in opposite clusters
%}

if (mod(n,2) == 0 ) % n is EVEN
%     generate Permutation Matrix, T
    I = eye(n);
    ix = randperm (n);
    T = I(ix,:);

%     generate adjaceny matrix A of a planted partition over n nodes
    n2 = n/2;
    P = random('bino', 1, p, n2, n2); % upper left block
    dP2 = random('bino', 1, p, n2, 1); % diagonal of the lower right block
    Q = random('bino', 1, q, n2, n2); % upper right block
%     carve the two triangular and diagonal matrices that we need
    U = triu(P, 1);
    L = tril(P,-1);
    dP = diag(P);
    B0 = U + U' + diag(dP);
    B1 = Q;
    B2 = Q';
    B3 = L + L' + diag(dP2);
    B =[B0 B1;B2 B3];

    comm1 = ones(1,n2);
    originalCluster = [comm1 -1*comm1];
    
% % % -----------------------------------------    
% %     % No RANDOMIZATION into CLUSTERS!
% %     A = B;
% %     partitionIndicatorVec  = originalCluster;
% % % % -----------------------------------------
%  PERMUTE THE NODES
%     Re-index Nodes of the graph. 
%     B*T' -> exchg columns; T*M -> exchg rows
    A = T*B*T'; 

%   Obtain the True_Cluster_NodeID for Graph A: Permute them
    partitionIndicatorVec  = originalCluster(ix);
% % -----------------------------------------
else
   warning('n == ODD!') ;
end
end





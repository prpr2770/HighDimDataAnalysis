% hw3q18: Not implemented!
clear all; close all; clc;

load zachary.mat

% From the question/visual graph: Identify true_Partition
nodeIDs = 1:34;
idx_teamA = [25 26 28 32 24 29 30 27 10 34 9 21 33 31 19 23 15 16]
idx_teamB = setdiff(1:34,idx_teamA)

truePartition = ones(size(nodeIDs));
truePartition(idx_teamB) = -1*ones(size(idx_teamB));


% Implement the algorithm to detect the partitions.

M = (A ~= 0 );
adjMatrix = zeros(size(A));
adjMatrix(M) = 1;

estPartition = runPartitionAlgo(adjMatrix);
cluster1 = (estPartition <0);
cluster1_idx_est = cluster1.*nodeIDs;

cluster2 = (estPartition >0);
cluster2_idx_est = cluster2.*nodeIDs;

overlapScore = getPartitionOverlap(truePartition, estPartition)





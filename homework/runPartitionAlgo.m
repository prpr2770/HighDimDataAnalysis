function partitionIndicatorVec = runPartitionAlgo(A)
%{ 
Q14.
Algorithm Partition
* compute the second dominant eigenvector, v2, of A, associated with the second largest
eigenvalue ?2.
* for i = 1 to n
if the coordinate i of v2 is positive,(v2)_i > 0, then  %what if its' ZERO?
assign node wi to community 1
else
assign node i to community 2.
end
end

partitionIndicatorVec : = {1(partition A), -1(partition B)}

%}


[V, D] = eigs(A,2);

vec = V(:,2)';
pos = (vec >0);
neg = (vec <0);


partitionIndicatorVec = pos - neg;



end

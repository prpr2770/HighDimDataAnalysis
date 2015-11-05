function  overlap = getPartitionOverlap(w1, w2)
%{
Q15. 
Derive the Overlap Metric based on the 
w1 : true partition vector
w2 : estimated partition vector
%}
n = length(w1);
del_w1_w2 = sum(w1 == w2);
del_minus_w1_w2 = sum(-w1 == w2);

rawoverlap = max(del_w1_w2, del_minus_w1_w2);

% random choice of w2 generates a non-zero overlap. Accounting for this:
overlap = (2/n)*rawoverlap - 1;         % Interpret: Prob. of successfully detecting communities.
end
 


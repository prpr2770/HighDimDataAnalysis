function [dist_metric dist_matrix] = runITMLonDataSet(X,y)
%{
Input:
X -> Data-vector coordinates: ROW-VEC
y -> SongGenres for each song: ROW-VEC
%}

dist_metric = MetricLearningAutotuneKnn(@ItmlAlg, y, X);
imagesc(dist_metric)
colorbar


% ------------------------------------------------------
% compute Symmetric Distance_Matrix
dist_matrix = getDistanceMatrixFromMetric(X,dist_metric);

end
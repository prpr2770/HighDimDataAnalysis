function DistanceMatrix = getDistanceMatrixFromMetric(X,distance_metric);

W = distance_metric;
totalSongs = length(X);
D = zeros(X,X);
for i = 1:totalSongs
    for j = i+1:totalSongs
    x_song = X(i,:);
    y_song = X(j,:);
    d_ij = (x_song - y_song)*W*(x_song - y_song)';
    if (d_ij >= 0)
        D(i,j) = d_ij.^(1/2);
    else
        warning('distance are complex')
    end
    
    end
end
Distance_Matrix = D + D';

end
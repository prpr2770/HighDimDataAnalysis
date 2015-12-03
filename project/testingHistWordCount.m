
X = rand(256,3);
tau = 10;
numCodeWords = length(X);
iterMax = 10;

hist_mat = zeros(numCodeWords,iterMax);
for i=1:iterMax
    
    numFrames = ceil(100*rand());
    Y = rand(numFrames,3);
    nbrs_of_songFrames = knnsearch(X,Y,'k',tau,'distance','euclidean');
    nbrs_of_songFrames = reshape(nbrs_of_songFrames,[],1);
    
    % extract histogram of occurence and re-structure into COL_VEC
    hist_vec = histc(nbrs_of_songFrames,1:numCodeWords);
    hist_vec = reshape(hist_vec,[],1);
    hist_mat(:,i) = (1/tau)*(1/numFrames)*hist_vec;
end

plot(hist_mat)
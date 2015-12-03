%{
Script to call all sequence of codes, to generate CodeWord Histograms and
TSNE Plots!

%}

% Extract Code-words from Batches of Songs
numCodeWordsPerBatch = 512; iter = 20; songsPerBatch = 10;
extract_MFCC_CodeWords_BatchOfSongs(numCodeWordsPerBatch,songsPerBatch,iter);

% Extract CodeWords for all the Batches-together. || ALLSONGS-CODEWORDS
numCodeWords = 512; tau = 10;iter = 50;
extract_MFCC_CodeWords_AllSongs(numCodeWords, iter)

% Compute CODEWORD-HIST for each song.
extract_MFCC_codeWord_Histograms( numCodeWords, tau);

% Execute DimensionalReduction
runCWhist_mfcc_tsne
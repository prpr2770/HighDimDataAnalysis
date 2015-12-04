%{
Script to do the following:

+> Read the tracks\ folder for all the songs. And generate the MFCC
information. Store in separate folder. 

+> extract the DYN_MFCC, FD_DYN_MFCC info from the MFCC extracted for each
song. 

+> Normalize the MFCC, DYN_MFCC and FD_MFCC and store info in song-wise
data.
Done by computing MEAN_[], STD_[] for each of above sets.( store this in
the aggregateData File! )

+>  Extract CodeWords from the given choice of data.  

% =========================================================================
Assumes songGenres information is contained in 
tracksDirName... : '\g1c_SongFeatures\songGenres.mat'

% data directory to store mfcc and dyn_mfcc of each song
matDataDirName = fullfile(tracksDirName,'song_data')

% data directory to store global/aggregate information
aggregateDataDirName = fullfile(tracksDirName,'aggregate_song_data')


% -------------------------------------------------------------------------------

% matfile to store aggData_allSongs
aggData_allSongs_fileName = strcat(aggregateDataDirName,'\aggData_allSongs.mat');
aggData_mFile = matfile(aggData_allSongs_fileName,'Writable',true);



%}

tracksDirName = 'H:\HighDimData\Project\ecen5322\Volumes\project\tracks\';

% -----------------------------------------------------------------------
% Extract Features
% Read .wav files and generate MFCC, FD_MFCC, DYN_MFCC
% extractDynMFCC_song_wise(tracksDirName);
% extractDynMFCC_song_wise


% extract Normalized MFCC, FD_MFCC, DYN_MFCC
% compute_songwise_normalized_data(tracksDirName);
% compute_songwise_normalized_data


% -----------------------------------------------------------------------
% Extract Code-words for DB:

dataType = 'MFCC'; % 'DYN_MFCC' 'FD_DYN_MFCC'

%{
% Extract from Batches of Songs
numCodeWordsPerBatch = 512; iter = 20; songsPerBatch = 10;
extract_codewords_batchwise(tracksDirName, numCodeWordsPerBatch,songsPerBatch,iter, dataType);
%}


% Extract CodeWords for all the Batches-together. || ALLSONGS-CODEWORDS
numCodeWords = 512; tau = 10;iter = 50;
% extract_MFCC_CodeWords_AllSongs(numCodeWords, iter, dataType)
extract_codewords_forallsongs(tracksDirName,numCodeWords, iter)


% -----------------------------------------------------------------------
% Compute CODEWORD-HIST for each song.
extract_codeword_histogram_allSongs(tracksDirName, numCodeWords, tau, dataType);

% =========================================================================
% %}


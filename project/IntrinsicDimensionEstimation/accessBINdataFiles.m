% Script to access the binary files
fileName = './Moebius/Moebius_20.BIN';
fileID = fopen(fileName);

% read the data
M = fread(fileID);

fclose(fileID);
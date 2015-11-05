function wavWindows = extractWindowsFromSong(wavFile, Fs)
%{
Input:
wavFile : row-vector of the audio-clip to be dissected
Fs : Sampling frequency

Output:
All the windows as rows of a matrix.
%}


% Obtain a song. And reshape the 1-D vector into a 2-D vector of 256 columns.
% Basically, generating overlapping intervals of songs - with each interval being 256 length.

% how to determine the cutting up of songs?
% 1. Determine the center data point.
% 2. first interval : where the center of song, is center of the 20ms interval concerned.
% 3. All other windows, grow from the center of the song, outwards.
% 4. the extremity windows, with datapoints less than the 20ms requirement, are rejected/thrown away!

% ------------------------------------------------------------------
% we want to cut songs into 100ms intervals? Fs = 11025
% 
% close all; clear all; clc
% 
% hfile2 = 'artist_1_album_1_track_1.wav';
% % hfile2 = 'artist_1_album_1_track_2.wav';
% % hfile2 = 'artist_78_album_1_track_3.wav';
% 
% wavInfo = audioinfo(hfile2)
% [wav, Fs] = audioread(hfile2);
% 
% wavFile = wav'; % Ensure ROW VECTOR
% size(wavFile)
% ------------------------------------------------------------------

intervalDuration = 0.04; %seconds = 200ms % Ensures that Odd number of points.
windowSize = ceil(intervalDuration*Fs); %221 samples!
halfWindow = (windowSize - 1)/2;
frameBreak = halfWindow + 1;

% 
% % ==============================================================
% % Code to test the Window-partitioning of WavFile
% wavFile = 1:1:11; % ensure ROW VECTOR!!!!
% 
% windowSize = 5; % always ODD!
% halfWindow = (windowSize-1)/2;
% frameBreak = halfWindow + 1;
% % ==============================================================

%find the central data-point:
lenWave = length(wavFile);
if mod(lenWave,2) ~= 0 % Odd length
    cntr = (lenWave + 1)/2;
else
    cntr = lenWave/2;
end

% split the signal at the center 
s_half = wavFile(cntr+1:end);
f_half = wavFile(1:cntr-1);

% ------------------------------------------------------
% examining second_half of song
numWindows = floor((length(s_half)+1)/frameBreak);
% reshape the s_half into a Matrix-of-frames
if length(s_half) >= (numWindows+1)*frameBreak
    s_half = s_half(1:(numWindows+1)*frameBreak );
elseif length(s_half) + 1 == numWindows*frameBreak
    s_half = [s_half 0];
else
    s_half = s_half(1:numWindows*frameBreak);
end

reShape = reshape(s_half,[frameBreak, numWindows]);
s_half_reshape = reShape';
clear reShape
M1 = s_half_reshape(1:end-1,:);
M2 = s_half_reshape(2:end,1:end-1);
s_half_Matrix = [M1 M2];

% ------------------------------------------------------
% examining central window
cntr_Matrix = wavFile(cntr - halfWindow : cntr + halfWindow);
% length(cntr_Matrix) == windowSize

% ------------------------------------------------------
% examining first_half of song

%mirror-image the first-half : Since the windows are computed
%outwards-from-center
f_half = fliplr(f_half);

%repeat method of structring
numWindows = floor((length(f_half)+1)/frameBreak);
% reshape the s_half into a Matrix-of-frames
if length(f_half) >= (numWindows+1)*frameBreak
    f_half = f_half(1:(numWindows+1)*frameBreak );
elseif length(f_half) + 1 == numWindows*frameBreak
    f_half = [f_half 0];
else
    f_half = f_half(1:numWindows*frameBreak);
end

reShape = reshape(f_half,[frameBreak, numWindows]);
f_half_reshape = reShape';
clear reShape
M1 = f_half_reshape(1:end-1,:);
M2 = f_half_reshape(2:end,1:end-1);
f_half_temp = [M1 M2];

f_half_Matrix = f_half_temp(end:-1:1,end:-1:1);

% ------------------------------------------------------
% Combining all parts to obtain - Windows
wavWindows = [f_half_Matrix; cntr_Matrix; s_half_Matrix];
% size(wavWindows)

end


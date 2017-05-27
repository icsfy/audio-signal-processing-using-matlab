%% An example of noise removal from an audio file

clearvars;
close all;

%% Read in the file
[f, Fs] = audioread('handel_J.wav');
T = 1 / Fs;             % Sampling period
L = length(f);          % Length of signal
t = (0:L-1) * T;        % Time vector

%% Plot audio channels
N = size(f, 1);
figure;  stem(t, f);
title('Original: Time-domain');  xlabel('time(seconds)');

%% Plot the spectrum
df = Fs / N;
w = (-(N/2):(N/2)-1)*df;
y = fft(f) / N;         % For normalizing
y2 = fftshift(y);
figure;  plot(w, abs(y2));
title('Original: Amplitude Spectrum');  xlabel('Frequency(Hz)');

%% plot the single-sided amplitude spectrum
figure;
plot(Fs*(0:(L/2))/L, abs(y2(N/2:end)))
title('Original: Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')

%% find f1
% [fmax1, f1index] = max(abs(y2(N/2:end)));
% fmax1, f1index = f1index*df

%% filtering f1: 1043Hz
n = 2;
beginFreq = 1040 / (Fs/2);
endFreq = 1050 / (Fs/2);
[b,a] = butter(n, [beginFreq, endFreq], 'stop');
fOut = filter(b, a, f);
% extract f1
[b,a] = butter(n, [beginFreq, endFreq], 'bandpass');
f1 = filter(b, a, f);
figure;  plot(w,abs(fftshift(fft(f1)/N)));  title('f1: 1043Hz');

%% find f2
% [fmax2, f2index] = max(abs(fftshift(fft(fOut)/N)));
% fmax2, f2index = Fs/2 - f2index*df

%% filtering f2: 1956Hz
n = 2;
beginFreq = 1950 / (Fs/2);
endFreq = 1960 / (Fs/2);
[b,a] = butter(n, [beginFreq, endFreq], 'stop');
fOut = filter(b, a, fOut);
% extract f2
[b,a] = butter(n, [beginFreq, endFreq], 'bandpass');
f2 = filter(b, a, f);
figure;  plot(w,abs(fftshift(fft(f2)/N)));  title('f2: 1956Hz');

%% For normalizing
Ym = max(max(max(fOut)),max(abs(min(fOut))));
fOut = fOut ./ Ym;

%% After processing
figure;  stem(t, fOut);
title('After processing: Time-domain');  xlabel('time(seconds)');
figure;  plot(w,abs(fftshift(fft(fOut)/N)));
title('After processing: Amplitude Spectrum');  xlabel('Frequency(Hz)');

%% Create object for playing audio
pOrig = audioplayer(f,Fs);  % original signal
% pOrig.play;
p = audioplayer(fOut, Fs);  % filtered signal
% p.play;
pf1 = audioplayer(f1, Fs);  % noise: f1: 1043Hz
% pf1.play;
pf2 = audioplayer(f2, Fs);  % noise: f2: 1956Hz
% pf2.play;

%% Write out to audio file
audiowrite('handel_J_processed.wav', fOut, Fs);
audiowrite('f1.wav',f1,Fs);
audiowrite('f2.wav',f2,Fs);


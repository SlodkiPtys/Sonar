%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: reflection.m
% Description: Matlab program to determine the reflection
% response and the homodyned wave
% Dependencies: Signal Processing Toolbox, Matlab v7.1 R14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The original wave sent out from the loudspeaker
original = audioread('original.wav');

% Convert the stereo 'original' to mono if it's stereo
if size(original, 2) == 2
    original = mean(original, 2);
end

% The reflection from the snowpack
recorded = audioread('recorded.wav');

% Convert the stereo 'recorded' to mono if it's stereo
if size(recorded, 2) == 2
    recorded = mean(recorded, 2);
end

% The length of the recorded signal
n = length(recorded);

% The air-coupled wave arriving on the side of the microphone
side = audioread('side.wav');

% The averaged environmental noise
en = audioread('noise.wav');

% Noise in the recording system found by the loopback test
nrc = audioread('nrc.wav');

% The sampling rate of the ADC (96000 Hz)
fs = 96000;
period = 1/fs;
t = period*(1:n);

% Filter all of the waves to prevent higher frequency components from being in the signal
% Load the filter coefficients (assuming b and a are the filter coefficients)
load('filter.mat'); % Load filter coefficients

% Filter the signals
original = filter(b, a, original);
recorded = filter(b, a, recorded);
side = filter(b, a, side);
en = filter(b, a, en);
nrc = filter(b, a, nrc);

% Assuming that the additional samples in 'recorded' are at the end and can be discarded:
recorded_trimmed = recorded(1:88200);

% Now perform the subtraction with the trimmed 'recorded' vector:
sPrime = recorded_trimmed - side - en - nrc;

% Assuming that the additional samples in 'original' are at the end and can be discarded:
original_trimmed = original(1:88200);

% Take trimmed original into the frequency domain
originalF = fft(original_trimmed);

% Take trimmed sPrime into the frequency domain
sPrimeF = fft(sPrime);

% Perform the division in the frequency domain
sDoublePrimeF = originalF./sPrimeF;

% Take the reflection response back into the time domain by use of the inverse FFT
sDoublePrimeT = ifft(sDoublePrimeF);
%Apply Weiner Spiking Deconvolution (WSD) to the reflection response
%to remove the snowpack attenuation
%set the length of the spiking operator
%this is set by trial and error
spikingop = 29;
%set the length of the pre-whitening
%this is also set by trial and error
prewhit = 0.1;
%size of the vector
n = size(sDoublePrimeT);
%compute the autocorrelation of the signals
ar = xcorr(sDoublePrimeT(:,1), sDoublePrimeT(:,1), spikingop);
%set up for the levinson recursion
s = ar(:,1).*hamming(spikingop*2+1);
ref = s(spikingop+1:2*spikingop,1);
p = ref(1,1)*(1+0.1/100.); %pre-whitening procedure
ref(1,1) = p;
%find the levinson coefficients
c = levinson(ref,spikingop);
c = c';
%convolution of the data with the filter
cdfilter = conv2(sDoublePrimeT,c);
dfg = cdfilter(1:n,:);
maxi = max(max(dfg));
mult = dfg * maxi/(max(max(sDoublePrimeT)));
%mult now contains the reflection response from the snowpack
%FMCW section of the code
%multiply both of the waves together (homodyning)
%note that what we do here is multiply the original wave
%the multiplication is repeated in the time domain
mm = sPrime .* original_trimmed;
%set up the parameters for the chirp-z transform
fs = 96000; f1 = 0.00001; f2 = 100; % in hertz
%fs = 96000; f1 = 20; f2 = 50;
m = 2048;
w = exp(-j*2*pi*(f2-f1)/(m*fs));
a = exp(j*2*pi*f1/fs);
%calculate the chirp-z transform
z = czt(mm,m,w,a);
N = length(z);
Pyy = z.* conj(z)/N;
Pyy = Pyy./1500000;
%The next step is to plot all of the data
% Plot the original wave sent out from the loudspeaker
subplot(2,2,1);
plot(original_trimmed);
title('Original Wave');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the recorded wave
subplot(2,2,2);
plot(recorded_trimmed);
title(' (b) Recorded Wave');
xlabel('Time (s)');
ylabel('Amplitude');
set(findobj(gca,'Type','line'), 'Color', 'black');

% Define time vector for the reflection coefficients (mult)
% This should match the length of 'mult' after any processing steps
t_mult = (0:length(mult)-1) / fs;

% Plot the reflection coefficients (mult)
subplot(2,2,3);
plot(t_mult, mult);
title(' (c) Reflection Coefficients');
xlabel('Time (s)');
ylabel('Magnitude');
set(findobj(gca,'Type','line'), 'Color', 'black');

% Plot the homodyned response in frequency domain
subplot(2,2,4);
plot(Pyy);
title(' (d) Homodyned Response in Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
set(findobj(gca,'Type','line'), 'Color', 'black');
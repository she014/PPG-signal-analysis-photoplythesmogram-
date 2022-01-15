% this matlab code includes several different PPG filtering approaches.
% the code and idea are demonstrated in "Analysis: an optimal filter for
% short photoplethysmogram signals." Liang et al., 2018
clc; clear all; close all;

% read a ppg signal, this signal was collected by a optical PPG sensor from
% fingertip, the sampling rate is 4000Hz.
ppg = csvread('ppg_sample.csv');
fs = 4000; % sampling rate
figure; plot(ppg); xlim([1 length(ppg)]); title('Raw PPG waveform'); xlabel('sample'); ylabel('amplitude');

% as it can be observed this PPG signal is noisy, includes powerline noise,
% baseline wander...

% step1: remove powerline noise, design and apply a 60Hz notch filter to
% remove powerline noise.
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
ppg_60 = filtfilt(d,ppg);
figure; hold all; plot(ppg,'b'); plot(ppg_60,'r'); xlim([1 length(ppg)]); title('Raw PPG waveform & removed powerline noise'); xlabel('sample');
legend('raw signal','removed 60Hz');
% it can be observed that there were not too much powerline noise, as two
% signals are basically the same. But sometimes there will be very strong
% powerline noises (50/60Hz). So it is necessary to apply a notch filter
% before further denoising.

fn = fs/2; % Shannon sampling rate which is 1/2 sampling frequency

% just to clarify. In the following examples, the "order" decides the
% property of the filter, I will use order=4 for all examples. "fl" is the
% low cutoff frequency, "fh" is high cutoff frequency. As for ppg, I will
% use fl = 0.5Hz, fh = 10Hz.

% method 1: butterworth bandpass filter
% [A,B,C,D] = butter(order,[fl fh]/fn);
[a1,b1,c1,d1] = butter(4,[0.5 10]/fn);
[filter_sos1,g1] = ss2sos(a1,b1,c1,d1);
ppg_butter = filtfilt(filter_sos1,g1,ppg_60);

% method 2: Chebyshev type1 bandpass filter
% [A,B,C,D] = cheby1(order,0.1,[fl fh]/fn)
[a2,b2,c2,d2] = cheby1(4,0.1,[0.5 10]/fn);
[filter_sos2,g2] = ss2sos(a2,b2,c2,d2);
ppg_ch1 = filtfilt(filter_sos2,g2,ppg_60);

% method 3: Chebyshev type2 bandpass filter
% [A,B,C,D] = cheby2(order,20,[fl fh]/fn)
[a3,b3,c3,d3] = cheby2(4,20,[0.5 10]/fn);
[filter_sos3,g3] = ss2sos(a3,b3,c3,d3);
ppg_ch2 = filtfilt(filter_sos3,g3,ppg_60);

% method 4: Elliptic filter (ellip)
% [A,B,C,D] = ellip(order,0.1,30,[fl fh]/fn)
[a4,b4,c4,d4] = ellip(4,0.1,30,[0.5 10]/fn);
[filter_sos4,g4] = ss2sos(a4,b4,c4,d4);
ppg_elp = filtfilt(filter_sos4,g4,ppg_60);

% method 5: fir hamming
% d = fir1(order,[fl fh]/fn,'bandpass');
d1 = fir1(4,[0.5 10]/fn,'bandpass');
ppg_hm = filtfilt(d1,1,ppg_60);

% method 6: fir least square method
% d = designfilt('bandpassfir','FilterOrder',order,'StopbandFrequency1',fL-0.2,'PassbandFrequency1',fL,'PassbandFrequency2',fH,'StopbandFrequency2',fH+2,'DesignMethod','ls','SampleRate',sample_freq);
d2 = designfilt('bandpassfir','FilterOrder',2,'StopbandFrequency1',0.5-0.2,'PassbandFrequency1',0.5,...
               'PassbandFrequency2',10,'StopbandFrequency2',10+2,'DesignMethod','ls','SampleRate',fs);
ppg_ls = filtfilt(d2,ppg_60);

% method 7: signal smoothing
ppg_sm = smooth(ppg_60, 200); % smooth the signal with a 200 samples window

% method 8: median filter
ppg_md = medfilt1(ppg_60, 200); % median filter with a window size of 200 samples, this can be adjusted 

% method 9: wavelet method
ppg_wav = wden(ppg_60,'modwtsqtwolog','s','mln',4,'db2');

% visualize the results
figure; set(gcf,'Position',[100 100 2000 800]);
subplot(3,3,1); hold all; plot(ppg_60,'b'); plot(ppg_butter,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','butterworth'); xlabel('samples'); ylabel('amplitude'); title('Butterworth bandpass');

subplot(3,3,2); hold all; plot(ppg_60,'b'); plot(ppg_ch1,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','chebyshev type1'); xlabel('samples'); ylabel('amplitude'); title('Chebyshev type1 bandpass');

subplot(3,3,3); hold all; plot(ppg_60,'b'); plot(ppg_ch2,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','chebyshev type2'); xlabel('samples'); ylabel('amplitude'); title('Chebyshev type2 bandpass');

subplot(3,3,4); hold all; plot(ppg_60,'b'); plot(ppg_elp,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','Elliptic'); xlabel('samples'); ylabel('amplitude'); title('Elliptic bandpass');

subplot(3,3,5); hold all; plot(ppg_60,'b'); plot(ppg_hm,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','FIR hamming'); xlabel('samples'); ylabel('amplitude'); title('FIR hamming bandpass');

subplot(3,3,6); hold all; plot(ppg_60,'b'); plot(ppg_ls,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','FIR least square'); xlabel('samples'); ylabel('amplitude'); title('FIR least square bandpass');

subplot(3,3,7); hold all; plot(ppg_60,'b'); plot(ppg_sm,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','smooth'); xlabel('samples'); ylabel('amplitude'); title('signal smoothing');

subplot(3,3,8); hold all; plot(ppg_60,'b'); plot(ppg_md,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','median filter'); xlabel('samples'); ylabel('amplitude'); title('Median filter');

subplot(3,3,9); hold all; plot(ppg_60,'b'); plot(ppg_wav,'r'); xlim([1 length(ppg_60)]); 
legend('before filter','wavelet'); xlabel('samples'); ylabel('amplitude'); title('wavelet');

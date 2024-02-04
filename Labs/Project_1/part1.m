clc;
clear;
close all;

del_t = 0.4 * 10^-3;        % this sampling period will be used to calculate the frequency range
n = linspace(-8, 7, 16);    % 9th sample is t=0
absTime = n * del_t;

sampled_time_v = [0.6875, 0.3261, 0.5429, 1.2022, ...
    1.6875, 1.6327, 1.3384, 1.3655, 1.8125, 2.1739, 1.9571, ...
    1.2978, 0.8125, 0.8673, 1.1616, 1.1345]; %9th sample is t=0 origin.

L = length(sampled_time_v);
fqs = 1/del_t * 1/L * (-L/2:L/2-1); % Specific sample / (sample period * Number of samples)
figure;
plot(absTime,sampled_time_v)    % plot time domain to visualize
title('Time Domain Data')
xlabel('Time (s)')
ylabel('Signal Strength')

fft = fft(fftshift(sampled_time_v)); %fftshift swaps left and right halves of vector. here, 9:16 get moved to front
figure;                              %thus array goes from 1:16 to 9:16then1:8   
subplot(2,1,1)
plot(fqs, fftshift(abs(fft)))                     %SOMEONE LABEL THESE PLOTS PLZ
subplot(2,1,2)
plot(fqs, fftshift(angle(fft)))% this should be the end of question 1 (with some explanation)
clear all;
close all;
del_t = 0.4 * 10^-3;
n = linspace(-8, 7, 16);

sampled_time_v = [0.6875, 0.3261, 0.5429, 1.2022, ...
    1.6875, 1.6327, 1.3384, 1.3655, 1.8125, 2.1739, 1.9571, ...
    1.2978, 0.8125, 0.8673, 1.1616, 1.1345]; %9th sample is origin. array was missing 1.6327
L = length(sampled_time_v);
fqs = 1/del_t * 1/L * (-L/2:L/2-1); %herp derp i looked it up. PLEASE CHANGE THIS
figure;
plot(n,sampled_time_v)

%Testing update speed
fft2 = fft(circshift(sampled_time_v,2));
fft1 = fft(sampled_time_v); %fftshift swaps left and right halves of vector. here, 9:16 get moved to front
figure;                              %thus array goes from 1:16 to 9:16then1:8   
subplot(2,1,1)
plot(fqs, abs(fft2), 'o')
hold on;
plot(fqs, abs(fft1));
hold on;
plot(fqs, fftshift(abs(fft1)))    
legend(" Circular shifted Time", "Base FFT", "FFTshifted")
%SOMEONE LABEL THESE PLOTS PLZ
subplot(2,1,2)
plot(fqs, angle(fft1))
hold on;
plot(fqs, fftshift(angle(fft1)))% this should be the end of question 1 (with some explanation)

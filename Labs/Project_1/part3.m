clc; clear; close all

%will generate two 1-D signals, then combine to a 2-D signal, then do 2D
%fft and two 1d fft and compare
%will generate a white gaussian and a sinusoid, each composed of 64 points from
baseValues = linspace(-2*pi, 2*pi, 128); % generating base signal from -2pi to 2pi. will use center -pi to pi portion
delta_t = baseValues(1)-baseValues(2);  % record sample rate

g_x = 1/sqrt(2*pi)*exp(-0.5*baseValues.^2); %generate & plot gaussian curve from -2pi to 2pi
subplot(8,4,1)
plot(baseValues,g_x)      
xlabel('t')
title('g_x')
subplot(8,4,3)
plot(baseValues(33:96), g_x(33:96)) %plot gaussian curve from -pi to pi
xlabel('t')
title('g_x base')

g_y = cos(baseValues);  %generate and plot cosine from -2pi to 2pi
subplot(8,4,2)
plot(baseValues,g_y)
xlabel('t')
title('g_y')
subplot(8,4,4)
plot(baseValues(33:96), g_y(33:96)) %plot cosine from -pi to pi
xlabel('t')
title('g_y base')

[tempX, tempY] = meshgrid(g_x, g_y);                    %generate fully sized array of X and Y to prep for 3D
[tempXX, tempYY] = meshgrid(g_x(33:96), g_y(33:96));    %same as above but limited to -pi to pi

g_f = tempX.*tempY;     % generate 2d plots for both limits and plot in 3d and in heat map for visualization
g_ff = tempXX.*tempYY;  % plotted both -2pi to 2pi and -pi to pi plots to imply the signals are repeating 
subplot(8,4,[5,6,9,10]) % and just a window of the repeating portion is being observed 
surf(baseValues, baseValues, g_f)
subplot(8,4,[13,14,17,18])
imagesc(baseValues, baseValues, g_f)
subplot(8,4,[7,8,11,12])
surf(baseValues(33:96), baseValues(33:96), g_ff)
subplot(8,4,[15,16,19,20])
imagesc(baseValues(33:96), baseValues(33:96), g_ff)

G_x = fft((g_x(33:96)));    % generate the 1D ffts using fftshift to maintin zero
G_y = fft((g_y(33:96)));    % position from the center of the array
subplot(8,4,21)
plot(fftshift(abs(G_x)))           % place ifftshift on frequency domain to show negative to positive freq 
subplot(8,4,22)                     % (mostly to show gaussian in, gaussian out)
plot(fftshift(abs(G_y)))
subplot(8,4,23)
plot(fftshift(angle(G_x)))
subplot(8,4,24)
plot(fftshift(angle(G_y)))

[fXmag,fYmag] = meshgrid(fftshift(abs(G_x)), fftshift(abs(G_y))); 
G_f_mag = fXmag.*fYmag;     %repeat process from start of problem to generate 2D fft as the 1D fft's multiplied
[fXang,fYang] = meshgrid(fftshift(angle(G_x)), fftshift(angle(G_y)));
G_f_angle = fXang.*fYang;
subplot(8,4,[25,26])
surf(G_f_mag)
subplot(8,4,[27,28])
surf(G_f_angle)

G_ff = fft2(g_ff);      %generate 2D fft directly from the 2D time domain plot for comparison
subplot(8,4,[29,30])
surf(fftshift(abs(G_ff)))
title('2-D fft Mag')
subplot(8,4,[31,32])
surf(fftshift(angle(G_ff)))
title('2-D fft ang')
% EE233 Lab 4 Matlab code
% 
%
% This program is used for calculate center frequcies in EE233 Lab 4
% All parameters (R1, R2, R3, R4, R5, gama, minfreq, maxfreq) could be set
% mannually in this program. One of C1 and C2 could be a vector for
% sweeping.
%
% Paramters:
% R1, R2, R3, R4, R5 and gama are the parameters corresponding to Figure 2
% in Lab 4
% minfreq and maxfreq are the minimum and maximum frequency in frequency
% domain analysis
%
% Units:
% resistor: ohm
% capacitor: Farad
% frequency: Hz
clear; clc; clf; close all;
% Set parameters
R1 = 240e3; % R1 = 240k ohm
R2 = 240e3; % R2 = 240k ohm
R3 = 2.4e3; % R3 = 2.4k ohm
R4 = 2.4e3; % R4 = 2.4k ohm
R5 = 100e3; % R5 = 100k ohm
gama = 0.25; % gama = 0.25
minfreq = 10; % frequency starts at 10 Hz
maxfreq = 1e6; % frequency ends at 5 kHz
% Set C1 and C2
C1 = 0.1e-6; % set C1 here
% C2 = linspace(0.001e-6,0.1e-6,1e2); % set C2 here
C2 = 0.01e-6;
if(length(C1) == 1 && length(C2) == 1) % C1 and C2 are both numbers
[f0, G0, fc1, fc2] = Lab4_filter_Gain(R1, R2, R3, R4, R5, C1, C2,...
gama, minfreq, maxfreq, 'SingleCalculation');
fprintf('Center frequency f0 = %6g Hz\n', f0);
fprintf('Gain at f0 is G0 = %6g \n', G0);
fprintf('Cutoff frequencies: fc1 = %6g Hz, fc2 = %6g Hz\n', fc1, fc2);
elseif (length(C1) > 1 && length(C2) == 1) % C1 is a vector
f0 = zeros(1, length(C1));
G0 = zeros(1, length(C1));
for m = 1: length(C1)
[f0(m), G0(m), fc1, fc2] ...
= Lab4_filter_Gain(R1, R2, R3, R4, R5, C1(m), C2,...
gama, minfreq, maxfreq, 'VectorCalculation');
end
plot(C1, f0);
xlabel('C1 (Farad)');
ylabel('Center frequency (Hz)');
grid on;
hold on;
elseif (length(C1) == 1 && length(C2) > 1) % C2 is a vector
f0 = zeros(1, length(C2));
G0 = zeros(1, length(C2));
for m = 1: length(C2)
[f0(m), G0(m), fc1, fc2] ...
= Lab4_filter_Gain(R1, R2, R3, R4, R5, C1, C2(m),...
gama, minfreq, maxfreq, 'VectorCalculation');
end
plot(C2, f0);
xlabel('C2 (Farad)');
ylabel('Center frequency (Hz)');
grid on;
hold on;
else
err = MException('Input:TooMuchInput',...
'Too much input.\nOnly one parameters could be swept.');
throw(err);
end
% Matlab sub-function for EE233 Lab 4 Filters
%
% This function is used to find
% center frequency, f0
% gain at center frequency, G0
% cutoff frequencies, fc1, fc2
% plot gain in frequency domain
% Given by
% Resistance (R1,R2,R3,R4,R5)
% Capacitance (C1,C5)
% Left portion of potentiometer (gama)
% frequency range
% Command: SingleCalculation, VectorCalculation
%
% IF R1 = R2 and R3 = R4
% If gama > 0.5, the gain is always less than 1, so it is a band-stop
filter
% If gama < 0.5, the gain is always greater than 1, so it is a band-pass
filter
% If gama = 0.5, the gain is always equal to 1
%
% Notice: Input references should all be considered as numbers
function [fc, Gc, fc1, fc2] = Lab4_filter_Gain(R1, R2, R3, R4, R5, C1, ...
C2, gama, minfreq, maxfreq, command)
% angular frequency range
samples = 1e6;
freq = linspace(minfreq, maxfreq, samples);
w = 2 * pi * freq;
% s domain in steady state
s = 1i * w;
% delta-Y transformation
Za = gama .* R5 ./ (1 + s * R5 * C1);
Zb = (1 - gama) * R5 ./ (1 + s * R5 * C1);
Zc = gama * (1 - gama) * R5 ^2 * C1 * s ./ (1 + s * R5 * C1);
% transfer function
Hs = (1 ./ (R3 + Za) + (Zc + 1 ./ (s .* C2)) ./ R1 .* ...
(1 ./ (R3 + Za) + 1./ (R4 + Zb) + 1./ (Zc + 1 ./ (s .* C2)))) ./ ...
(1 ./ (R4 + Zb) + (Zc + 1 ./ (s .* C2)) ./ R2 .* ...
(1 ./ (R3 + Za) + 1./ (R4 + Zb) + 1./ (Zc + 1 ./ (s .* C2))));
% Gain and its plot
gain = abs(Hs);
if(strcmp(command, 'SingleCalculation'))
clf;
loglog(freq, gain);
hold on;
grid on;
xlabel('frequency (Hz)');
ylabel('gain');
end
% find the range of gain
max_gain = max(gain);
min_gain = min(gain);

% find the center frequency and its gain
if (gama > 0.5) % band-stop filter
fc = freq(find(gain == min_gain, 1, 'first'));
Gc = min_gain;
elseif (gama < 0.5) % band-pass filter
fc = freq(find(gain == max_gain, 1, 'first'));
Gc = max_gain;
else % Hs = 1
err = MException ('Gain:ConstantGain', ...
'The gain is always 1 and the center frequency cannot be
identified. Please set gama not equal to 0.5');
throw(err);
end
% find the cutoff frequencies
if(strcmp(command, 'SingleCalculation'))
if ((gama > 0.5 && max_gain >= (1 / sqrt(2))) || (gama < 0.5 &&
max_gain <= sqrt(2)))
errBound = MException('CutoffFrequency:NoCutoffFrequency',...
'The gain is too small that there is no cutoff frequecies for
this filter. Please change to other options.');
throw(errBound);
else
cutoff_gain = max_gain / sqrt(2);
if (gama > 0.5) % band-stop filter
fc1 = freq(find(gain <= cutoff_gain, 1, 'first'));
fc2 = freq(find(gain <= cutoff_gain, 1, 'last'));
else % band-pass filter
fc1 = freq(find(gain >= cutoff_gain, 1, 'first'));
fc2 = freq(find(gain >= cutoff_gain, 1, 'last'));
end
if (isempty(fc1) || isempty(fc2)... % the frequency range
is not enough to cover the cutoff frequency
|| fc1 == fc2 || fc1 == 0 || fc2 == maxfreq)
errBound = MException('CutoffFrequency:OutOfBound',...
'The frequency range is not enough to cover two cutoff
frequencies. Please increase your frequency range');
throw(errBound);
end
end
else
fc1 = 0;
fc2 = 0;
end
% title on the plot
str = sprintf('Center frequency f0 = %2f Hz\n Gain at f0 is G0 = %2f \n
Cutoff frequencies fc1 = %2f Hz, fc2 = %2f Hz', ...
fc, Gc, fc1, fc2);
title(str);
end

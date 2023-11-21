clc
clear all
close all

%% Signal generation
fs=10e3;                                   % sampling frequency
A=1;                                       % amplitude
f1=1540;                                   % starting frequency
f2=1920;                                   % ending frequency
tau=2;                                     % signal duration
k0=(f2-f1)/(tau*f1*f2);                    % period slope
t=(0:1/fs:(tau-1/fs));
xn=A*exp(-1i*2*pi/k0*log(-k0*t+1/f1));     % pure signal
N=length(xn);                              % the number of discrete samples N (i.e., length)

%% Noise generation
snr=-5;                                    % SNR defined in decibels as 10log10(A2/σ2).
sigma=A*sqrt(1/(10^(snr/10)));             % variance of the noise
rng(2023);                                 % fixed noise seed
wn=sigma*randn(1,N);                       % white Gaussian noise
zn=xn+1.0*wn;                              % noisy signal

%% Parameters estimation
[k0estirefinement,f1estirefinement]=GD_IRLS(zn,fs);
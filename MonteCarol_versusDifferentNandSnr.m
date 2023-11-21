clc
clear all
close all

LINE_WIDTH=1.5;

%% Signal generation
fs=10e3;                           % sampling frequency
A=1;                               % amplitude
f1=1540;                           % starting frequency
f2=1920;                           % ending frequency
tauall=[0.25,0.5,1,2];             % signal duration (corresponding to the number of discrete samples N = 2500, 5000, 10000, and 20000)
Ntau=length(tauall);
k0all=(f2-f1)./(tauall*f1*f2);     % period slope (corresponding to the number of discrete samples N = 2500, 5000, 10000, and 20000)

% the pure signal 1 with the number of discrete samples N = 2500
tau1 = tauall(1);
t1 = (0:1/fs:(tau1-1/fs));
xn1 = A*exp(-1i*2*pi/k0all(1)*log(-k0all(1)*t1+1/f1));
N1 = length(xn1);

% the pure signal 2 with the number of discrete samples N = 5000
tau2 = tauall(2);
t2 = (0:1/fs:(tau2-1/fs));
xn2 = A*exp(-1i*2*pi/k0all(2)*log(-k0all(2)*t2+1/f1));
N2 = length(xn2);

% the pure signal 3 with the number of discrete samples N = 10000
tau3 = tauall(3);
t3 = (0:1/fs:(tau3-1/fs));
xn3 = A*exp(-1i*2*pi/k0all(3)*log(-k0all(3)*t3+1/f1));
N3 = length(xn3);

% the pure signal 4 with the number of discrete samples N = 20000
tau4 = tauall(4);
t4 = (0:1/fs:(tau4-1/fs));
xn4 = A*exp(-1i*2*pi/k0all(4)*log(-k0all(4)*t4+1/f1));
N4 = length(xn4);

%% Monte Carol simulation 
snrall=-10:1:0;          % SNR defined in decibels as 10log10(A2/σ2)
NSNR=length(snrall);
Np=2000;                 % 2000 Monte Carlo simulation runs

f1eall_refinement=zeros(Ntau,NSNR,Np);
k0eall_refinement=zeros(Ntau,NSNR,Np);

for kk=1:NSNR
    snr=snrall(kk);
    [snr, kk, NSNR]
    sigma=A*sqrt(1/(10^(snr/10)));
    for ll=1:Np
        % Noise generation
        % N = 2500
        wn1 = sigma*randn(1,N1);
        zn1 = xn1+1.0*wn1;

        % N = 5000
        wn2 = sigma*randn(1,N2);
        zn2 = xn2+1.0*wn2;

        % N = 10000
        wn3 = sigma*randn(1,N3);
        zn3 = xn3+1.0*wn3;

        % N = 20000
        wn4 = sigma*randn(1,N4);
        zn4 = xn4+1.0*wn4;

        % Parameters estimation and results saving
        % N = 2500
        [k0estirefinement1,f1estirefinement1] = GD_IRLS(zn1,fs);
        f1eall_refinement(1,kk,ll)=f1estirefinement1;
        k0eall_refinement(1,kk,ll)=k0estirefinement1;
        
        % N = 5000
        [k0estirefinement2,f1estirefinement2] = GD_IRLS(zn2,fs);
        f1eall_refinement(2,kk,ll)=f1estirefinement2;
        k0eall_refinement(2,kk,ll)=k0estirefinement2;

        % N = 10000
        [k0estirefinement3,f1estirefinement3] = GD_IRLS(zn3,fs);
        f1eall_refinement(3,kk,ll)=f1estirefinement3;
        k0eall_refinement(3,kk,ll)=k0estirefinement3;

        % N = 20000
        [k0estirefinement4,f1estirefinement4] = GD_IRLS(zn4,fs);
        f1eall_refinement(4,kk,ll)=f1estirefinement4;
        k0eall_refinement(4,kk,ll)=k0estirefinement4;
    end
end

%% Normalized root mean square error (NRMSE)
% the starting frequency f1
f1se_refinement=(f1eall_refinement-f1).*(f1eall_refinement-f1);
f1nrmse_refinement=sqrt(mean(f1se_refinement,3))/f1;

figure(1)
plot(snrall,log10(f1nrmse_refinement(1,:)),'>-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(f1nrmse_refinement(2,:)),'o-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(f1nrmse_refinement(3,:)),'x-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(f1nrmse_refinement(4,:)),'s-','LineWidth',LINE_WIDTH)
hold off
legend('{\itN}=2500','{\itN}=5000','{\itN}=10000','{\itN}=20000',fontsize=18,FontName='Times New Roman')
legend('boxoff')
grid
ylabel('Log10({\itNRMSE}) of {\itf}_1 (dB)',fontsize=18,FontName='Times New Roman')
set(gca,fontsize=18,FontName='Times New Roman')

% the period slope k0
k0se_refinement=zeros(Ntau,NSNR,Np);
k0nrmse_refinement=zeros(Ntau,NSNR);
for k=1:Ntau
    k0se_refinement(k,:,:)=(k0eall_refinement(k,:,:)-k0all(k)).*(k0eall_refinement(k,:,:)-k0all(k));
    k0nrmse_refinement(k,:)=sqrt(mean(k0se_refinement(k,:,:),3))/abs(k0all(k));
end

figure(2)
plot(snrall,log10(k0nrmse_refinement(1,:)),'>-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(k0nrmse_refinement(2,:)),'o-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(k0nrmse_refinement(3,:)),'x-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(k0nrmse_refinement(4,:)),'s-','LineWidth',LINE_WIDTH)
hold off
legend('{\itN}=2500','{\itN}=5000','{\itN}=10000','{\itN}=20000',fontsize=18,FontName='Times New Roman')
legend('boxoff')
grid
ylabel('Log10({\itNRMSE}) of {\itk}_0 (dB)',fontsize=18,FontName='Times New Roman')
set(gca,fontsize=18,FontName='Times New Roman')

%% Normalized mean absolute error (NMAE)
% the starting frequency f1
f1ae_refinement=abs(f1eall_refinement-f1);
f1nmae_refinement=(mean(f1ae_refinement,3))/f1;

figure(3)
plot(snrall,log10(f1nmae_refinement(1,:)),'>-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(f1nmae_refinement(2,:)),'o-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(f1nmae_refinement(3,:)),'x-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(f1nmae_refinement(4,:)),'s-','LineWidth',LINE_WIDTH)
hold off
legend('{\itN}=2500','{\itN}=5000','{\itN}=10000','{\itN}=20000',fontsize=18,FontName='Times New Roman')
legend('boxoff')
grid
ylabel('Log10({\itNMAE}) of {\itf}_1 (dB)',fontsize=18,FontName='Times New Roman')
set(gca,fontsize=18,FontName='Times New Roman')

% the period slope k0
k0ae_refinement=zeros(Ntau,NSNR,Np);
k0nmae_refinement=zeros(Ntau,NSNR);
for k=1:Ntau
    k0ae_refinement(k,:,:)=abs(k0eall_refinement(k,:,:)-k0all(k));
    k0nmae_refinement(k,:)=(mean(k0ae_refinement(k,:,:),3))/abs(k0all(k));
end

figure(4)
plot(snrall,log10(k0nmae_refinement(1,:)),'>-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(k0nmae_refinement(2,:)),'o-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(k0nmae_refinement(3,:)),'x-','LineWidth',LINE_WIDTH)
hold on
plot(snrall,log10(k0nmae_refinement(4,:)),'s-','LineWidth',LINE_WIDTH)
hold off
legend('{\itN}=2500','{\itN}=5000','{\itN}=10000','{\itN}=20000',fontsize=18,FontName='Times New Roman')
legend('boxoff')
grid
ylabel('Log10({\itNMAE}) of {\itk}_0 (dB)',fontsize=18,FontName='Times New Roman')
set(gca,fontsize=18,FontName='Times New Roman')

return
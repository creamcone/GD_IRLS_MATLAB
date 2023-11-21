function [k0esti,f1esti]=GD_IRLS(zn,fs)
%% The description of the algorithm
% (1) The function
%     By the algorithm, the frequency parameters of HFM signal, including the starting frequency and period slope, can be estimated based on the characteristics of group delay (GD).
%     The bandwidth of the HFM signal is preliminarily estimated by using the histogram statistics, 
%     and it is accurately estimated by using the transient features of amplitude spectrum,
%     the result of it is called the frequency distribution range (FDR).
%     The GD in the FDR of the HFM signal is extracted,
%     and the frequency parameters are estimated by combining the continuity criterion and iteratively reweighted least squares (IRLS) linear fitting on the GD.
% (2) The input
%     zn:      the analyzed HFM signal
%     fs:      sampling frequency
% (3) The output
%     k0esti:  the estimated period slope
%     f1esti:  the estimated starting frequency

%% STEP 1 Estimate the frequency distribution range (FDR) and extract the group delay (GD) of the HFM signal
%% STEP 1 A
%  Obtain the amplitude spectrum (AS) & spectrum phase (SP)
N=length(zn);
deltaf=fs/N;                 % frequency resolution
xf=fft(zn);
xfa=abs(xf(1:N/2));          % AS
xfp=angle(xf(1:N/2));        % SP

%% STEP 1 B
%  Normalize and iteratively smooth the AS
LF=length(xfa);
fdw=fs/20;
indexdw=round(fdw/deltaf);
xfa([1:indexdw LF-indexdw+1:LF])=0;

xfas0=xfa;
sumold=0;
for k=1:10
    xfasn=smooth(xfas0,9);               % smoothed AS (SAS)
    sumnew=sum(abs(xfasn'-xfa).^2);
    if  k>1 && abs(sumnew-sumold)<0.01*sumnew
        break;
    end
    sumold=sum(abs(xfasn'-xfa).^2);
    xfas0=xfasn;
end
xfais=xfasn;                             % iterative SAS (ISAS)
xfanis=xfais/max(xfais);                 % normalized ISAS (NISAS)
xfanis=xfanis';        

%% STEP 1 C
%  Perform the histogram statistics
[~,indexsM]=max(xfanis);
xc=0.05:0.1:0.95;
[H,~]=hist(xfanis,xc);

%% STEP 1 D
%  Estimate the bandwidth
H=H/max(H);
index=find((xc<0.4)|(xc>0.9));                 % the considered high-energy statistical range is [0.4,0.9]
H(index)=0;
[~,mh]=max(H);
lambda=xc(mh);

indexsig=find(xfanis>max(lambda-0.1,0.7));     % preliminarily determine the signal bandwidth according to the adaptive threshold and minimum threshold (0.7 Corresponds to 3dB bandwidth)
indexsigs=sort(indexsig);
nh=length(indexsigs);                          % bandwidth

%% STEP 1 E
%  Estimate the gravity frequency (GF)
indexsigsel=indexsigs([min(3,nh):max(nh-3,1)]);
if length(indexsigsel)<3
   indexsigsel=indexsigs;
end
kg=sum(indexsigsel.*xfanis(indexsigsel))/sum(xfanis(indexsigsel));
kg=round(kg);
GF=(kg-1)*deltaf;                     % GF

%% STEP 1 F
%  Calculate the upward feature (UF) and downward feature (DF)
sigmak=diff(xfanis);
sigmak=[sigmak,0];

krwmul=max(round(0.2*nh),round(0.05*kg));
if rem(krwmul,2)==0
    krwmul=krwmul+1;
end
krw=max((krwmul-1)/2,1);                         % the reference window length

Crk=zeros(1,LF);
Cfk=zeros(1,LF);
for k=krw+1:LF-1-krw
    kd=k-krw:k+krw;
    sigmakkd=sigmak(kd);
    Crk(k)=length(find(sigmakkd>0));
    Cfk(k)=length(find(sigmakkd<0));
end
Crk=Crk/max(Crk);
Cfk=Cfk/max(Cfk);

Uk=sigmak.*Crk;
Dk=sigmak.*Cfk;
Uk=Uk/max(abs(Uk));        % UF
Dk=Dk/max(abs(Dk));        % DF

%% STEP 1 G
%  Estimate the lower bound frequency (LBF) and upper bound frequency (UBF)
ksw=round(1.1*min(round(kg/5),nh));         % the searching window length

[~,indexLL]=max(Uk(max(kg-1,1):-1:max(kg-ksw,1)));
klb=kg-indexLL;
[~,indexRR]=min(Dk(min(kg+1,LF):1:min(kg+ksw,LF)));
khb=kg+indexRR;
if isempty(klb)==1
    klb=kg;
end
if isempty(khb)==1
    khb=kg;
end

flb=(klb-1)*deltaf;                       % LBF
fhb=(khb-1)*deltaf;                       % UBF
fpeak=(indexsM-1)*deltaf;                 % the peak frequency

%% STEP 1 H
%  Extract the GD in the FDR
xfpd=diff(xfp);
for k=1:N/2-1
    kk=0;
    while xfpd(k)>0 && kk<10
        kk=kk+1;
        xfpd(k)= xfpd(k)-2*pi;  
    end
end
taudka=-xfpd((1:N/2-1))/(2*pi*deltaf);           % GD in the full band

fk=((klb:khb)-1)*deltaf;
taudk=taudka(klb:khb);                           % GD in the FDR
xk=1./fk;
yk=taudk;

%% STEP 2 Eliminate the outliers and perform the IRLS linear fitting on the extracted GD
%% STEP 2 A
%  To eliminate the isolated outliers by the criterion of continuity
xkm=xk;
ykm=yk;
LA=length(ykm);
IndexOutlier=[];                            % isolated outliers
maximumtaudk=max(ykm);

for i=3:LA-2
    if  (      (abs(ykm(i)-ykm(i-1))>0.25*maximumtaudk && abs(ykm(i)-ykm(i+1))>0.1*maximumtaudk)...
            || (abs(ykm(i)-ykm(i+1))>0.25*maximumtaudk && abs(ykm(i)-ykm(i-1))>0.1*maximumtaudk)...
            || (abs(ykm(i)-ykm(i-1))>0.15*maximumtaudk && abs(ykm(i)-ykm(i+1))>0.15*maximumtaudk)...
            || (abs(ykm(i)-ykm(i-2))>0.25*maximumtaudk && abs(ykm(i)-ykm(i+2))>0.1*maximumtaudk)...
            || (abs(ykm(i)-ykm(i+2))>0.25*maximumtaudk && abs(ykm(i)-ykm(i-2))>0.1*maximumtaudk)...
            || (abs(ykm(i)-ykm(i-2))>0.15*maximumtaudk && abs(ykm(i)-ykm(i+2))>0.15*maximumtaudk)...
            || (abs(ykm(i)-ykm(i-2))>0.25*maximumtaudk && abs(ykm(i)-ykm(i+1))>0.1*maximumtaudk)...
            || (abs(ykm(i)-ykm(i+2))>0.25*maximumtaudk && abs(ykm(i)-ykm(i-1))>0.1*maximumtaudk))
        IndexOutlier=[IndexOutlier i];
    end
end

xkmnormal=xkm;
ykmnormal=ykm;
xkmnormal(IndexOutlier)=[];
ykmnormal(IndexOutlier)=[];                      % normal GD points

%% STEP 2 B & C
%  Perform the IRLS linear fitting on the extracted GD
[b]=robustfit(xkmnormal,ykmnormal);

%% STEP 3 Estimate the frequency parameters
k0esti=-1/b(2);
f1esti=-b(2)/b(1);

% The estimated f1 (f1esti) cannot be negative
% Therefore, the peak frequency is taken if f1esti is negative
if f1esti<0
    f1esti=fpeak;
end

end

function [b]=robustfit(x,y)
%% The description of the IRLS algorithm
% (1) The function
%     By the algorithm, the iteratively reweighted least squares (IRLS) linear fitting can be performed.
%     The maximum number of iterations is 10 and 
%     when the relative error of two iterations is less than 0.01, the iteration is stopped.
%     For each iteration, the weight less than 0.75 is set to 0,
%     to further reduce the impact of consecutive outliers on parameter estimation.
% (2) The input
%     x:       the argument for fitting the line
%     y:       the dependent variable for fitting the line
% (3) The output
%     p:       the parameters of the fitted line

%%
x_fit=x';
y_fit=y';
n_fit=length(x_fit);
G=[ones(1,n_fit);x_fit'].';

% the conventional least squares (CLS)
b_0=(G'*G)\G'*y_fit;
ytfit_0=b_0(2)*x_fit+b_0(1);
JLast_0=sum(abs(y_fit-ytfit_0).^2);

% the iteratively reweighted linear fitting (IRLS)
Qmax=10;
q=1;
JLast=JLast_0;

while(1)
    dgt=abs(y_fit-ytfit_0);
    w=(max(dgt)-dgt)./(max(dgt)-min(dgt));
    index=find(w<0.75);
    w(index)=0;
    W_q=diag(w);
    b_q=(G'*W_q*G)\(G')*W_q*y_fit;
    ytfit_q=b_q(2)*x_fit+b_q(1);
    JCur=sum(w.*abs(y_fit-ytfit_q).^2);
    if (abs(JCur-JLast)/JLast<0.01) || (q==Qmax)
        break;
    end
    JLast=JCur;
    ytfit_0=ytfit_q;
    q=q+1;
end
b=b_q;
end
function [k0esti,f1esti]=GD_IRLS(xt,fs)
%% The description of the algorithm
% (1) The function
%     By the algorithm, the frequency parameters of HFM signal, including the starting frequency and period slope, can be estimated based on the characteristics of group delay (GD).
%     The bandwidth of the HFM signal is preliminarily estimated by using the histogram statistics, 
%     and it is accurately estimated by using the transient features of amplitude spectrum,
%     the result of it is called the frequency distribution range (FDR).
%     The GD in the FDR of the HFM signal is extracted,
%     and the frequency parameters are estimated by combining the continuity criterion and iteratively reweighted least squares (IRLS) linear fitting on the GD.
% (2) The input
%     xt:      the analyzed HFM signal
%     fs:      sampling frequency
% (3) The output
%     k0esti:  the estimated period slope
%     f1esti:  the estimated starting frequency

%% STEP 1 Estimate the frequency distribution range (FDR) and extract the group delay (GD) of the HFM signal
%% STEP 1 A
%  Obtain the amplitude spectrum (AS) & spectrum phase (SP)
L=length(xt);
deltaf=fs/L;                 % frequency resolution
fx=(0:L/2-1)*deltaf;
xtf=fft(xt);
xfa=abs(xtf(1:L/2));         % AS
xfp=angle(xtf(1:L/2));       % SP

%% STEP 1 B
%  Normalize and iteratively smooth the AS
LF=length(xfa);
fdw=fs/20;
indexdw=round(fdw/deltaf);
xfa([1:indexdw LF-indexdw+1:LF])=0;

XFRO=xfa;
sumold=0;
for k=1:10
    XFRN=smooth(XFRO,9);               % smoothed AS (SAS)
    sumnew=sum(abs(XFRN'-xfa).^2);
    if  k>1 && abs(sumnew-sumold)<0.01*sumnew
        break;
    end
    sumold=sum(abs(XFRN'-xfa).^2);
    XFRO=XFRN;
end
XFRS=XFRN;                             % iterative SAS (ISAS)
xfasn=XFRS/max(XFRS);                  % normalized ISAS (NISAS)
xfasn=xfasn';        

%% STEP 1 C
%  Perform the histogram statistics
[~,indexsM]=max(xfasn);
xx=0.05:0.1:0.95;
[nh,xh]=hist(xfasn,xx);

%% STEP 1 D
%  Estimate the bandwidth
nhn=nh/max(nh);
nhr=nhn;
index=find((xh<0.4)|(xh>0.9));                 % the considered high-energy statistical range is [0.4,0.9]
nhr(index)=0;
[~,mh]=max(nhr);
Amp_Sig=xx(mh);

indexsig=find(xfasn>max(Amp_Sig-0.1,0.7));     % preliminarily determine the signal bandwidth according to the signal amplitude and minimum threshold (0.7 Corresponds to 3dB bandwidth)
indexsigs=sort(indexsig);
FB_LS=length(indexsigs);                       % bandwidth

%% STEP 1 E
%  Estimate the gravity frequency (GF)
indexsigsel=indexsigs([min(3,FB_LS):max(FB_LS-3,1)]);
if length(indexsigsel)<3
   indexsigsel=indexsigs;
end
indexweight=sum(indexsigsel.*xfasn(indexsigsel))/sum(xfasn(indexsigsel));
indexweight=round(indexweight);
fw=(indexweight-1)*deltaf;                     % GF

%% STEP 1 F
%  Calculate the upward feature (UF) and downward feature (DF)
xfasnd=diff(xfasn);
xfasnd=[xfasnd,0];

wcl=max(round(0.2*FB_LS),round(indexweight*0.05));
if rem(wcl,2)==0
    wcl=wcl+1;
end
wclh=max((wcl-1)/2,1);                         % the reference window length

xfasnup=zeros(1,LF);
xfasnde=zeros(1,LF);
for kk=wclh+1:LF-1-wclh
    indexkk=kk-wclh:kk+wclh;
    tempkkd=xfasnd(indexkk);
    indexp=find(tempkkd>0);
    indexn=find(tempkkd<0);
    xfasnup(kk)=length(indexp);
    xfasnde(kk)=length(indexn);
end
xfasnup=xfasnup/max(xfasnup);
xfasnde=xfasnde/max(xfasnde);

xfasnupmd=xfasnd.*xfasnup;
xfasndemd=xfasnd.*xfasnde;
xfasnupmd=xfasnupmd/max(abs(xfasnupmd));        % UF
xfasndemd=xfasndemd/max(abs(xfasndemd));        % DF

%% STEP 1 G
%  Estimate the lower bound frequency (LBF) and upper bound frequency (UBF)
indexw=round(1.1*min(round(indexweight/5),FB_LS));         % the searching window length

[~,indexLL]=max(xfasnupmd(max(indexweight-1,1):-1:max(indexweight-indexw,1)));
indexsL=indexweight-indexLL;
[~,indexRR]=min(xfasndemd(min(indexweight+1,LF):1:min(indexweight+indexw,LF)));
indexsR=indexweight+indexRR;
if isempty(indexsL)==1
    indexsL=indexweight;
end
if isempty(indexsR)==1
    indexsR=indexweight;
end

fL=(indexsL-1)*deltaf;                 % LBF
fH=(indexsR-1)*deltaf;                 % UBF
fM=(indexsM-1)*deltaf;                 % the peak frequency

%% STEP 1 H
%  Extract the GD in the FDR
xfpdu=diff(xfp);
for k=1:L/2-1
    kk=0;
    while xfpdu(k)>0 && kk<10
        kk=kk+1;
        xfpdu(k)= xfpdu(k)-2*pi;  
    end
end
taudka=-xfpdu((1:L/2-1))/(2*pi*deltaf);           % GD in the full band

fk=((indexsL:indexsR)-1)*deltaf;
taudk=taudka(indexsL:indexsR);                    % GD in the FDR
xk=1./fk;
yk=taudk;

%% STEP 2 Eliminate the outliers and perform the IRLS linear fitting on the extracted GD
%% STEP 2 A
%  To eliminate the isolated outliers by the criterion of continuity
xkm=xk;
ykm=yk;
LTM=length(ykm);
IndexOutlier=[];                            % isolated outliers
maximum=max(ykm);

for mm=3:LTM-2
    if  (      (abs(ykm(mm)-ykm(mm-1))>0.25*maximum && abs(ykm(mm)-ykm(mm+1))>0.1*maximum)...
            || (abs(ykm(mm)-ykm(mm+1))>0.25*maximum && abs(ykm(mm)-ykm(mm-1))>0.1*maximum)...
            || (abs(ykm(mm)-ykm(mm-1))>0.15*maximum && abs(ykm(mm)-ykm(mm+1))>0.15*maximum)...
            || (abs(ykm(mm)-ykm(mm-2))>0.25*maximum && abs(ykm(mm)-ykm(mm+2))>0.1*maximum)...
            || (abs(ykm(mm)-ykm(mm+2))>0.25*maximum && abs(ykm(mm)-ykm(mm-2))>0.1*maximum)...
            || (abs(ykm(mm)-ykm(mm-2))>0.15*maximum && abs(ykm(mm)-ykm(mm+2))>0.15*maximum)...
            || (abs(ykm(mm)-ykm(mm-2))>0.25*maximum && abs(ykm(mm)-ykm(mm+1))>0.1*maximum)...
            || (abs(ykm(mm)-ykm(mm+2))>0.25*maximum && abs(ykm(mm)-ykm(mm-1))>0.1*maximum))
        IndexOutlier=[IndexOutlier mm];
    end
end

xkmcd=xkm;
ykmcd=ykm;
xkmcd(IndexOutlier)=[];
ykmcd(IndexOutlier)=[];                      % normal GD points

%% STEP 2 B & C
%  Perform the IRLS linear fitting on the extracted GD
[pr]=robustfit(xkmcd,ykmcd);

%% STEP 3 Estimate the frequency parameters
k0esti=-1/pr(2);
f1esti=-pr(2)/pr(1);

% The estimated f1 (f1esti) cannot be negative
% Therefore, the peak frequency (fM) is taken if f1esti is negative
if f1esti<0
    f1esti=fM;
end

end

function [p]=robustfit(x,y)
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
p_0=(G'*G)\G'*y_fit;
ytfit_0=p_0(2)*x_fit+p_0(1);
JLast_0=sum(abs(y_fit-ytfit_0).^2);

% the iteratively reweighted linear fitting (IRLS)
Nmax=10;
i=1;
JLast=JLast_0;

while(1)
    dgt=abs(y_fit-ytfit_0);
    w=(max(dgt)-dgt)./(max(dgt)-min(dgt));
    index=find(w<0.75);
    w(index)=0;
    W_q=diag(w);
    p_q=(G'*W_q*G)\(G')*W_q*y_fit;
    ytfit_q=p_q(2)*x_fit+p_q(1);
    JCur=sum(w.*abs(y_fit-ytfit_q).^2);
    if (abs(JCur-JLast)/JLast<0.01) || (i==Nmax)
        break;
    end
    JLast=JCur;
    ytfit_0=ytfit_q;
    i=i+1;
end
p=p_q;
end
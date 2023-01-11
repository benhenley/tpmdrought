% test_tpmdrought.m 
% Simple script to demonstrate tpmdrought.m
% Code by Ben Henley: @benhenley, www.benhenley.net
%
% Code history: 
% v1.00: Initial commit
%
% Disclaimer:
% This software is provided “as is”, without warranty of any kind, express 
% or implied, including but not limited to the warranties of merchantability, 
% fitness for a particular purpose and noninfringement. In no event shall 
% the authors or copyright holders be liable for any claim, damages or other 
% liability, whether in an action of contract, tort or otherwise, arising 
% from, out of or in connection with the software or the use or other dealings 
% in the software.
%

%% 1. Set up ---- 
clear;clc;close all;format compact;
set(0,'defaultAxesFontSize',18);set(0,'DefaultLegendAutoUpdate','off')
rng('default'); %for reproducibility

%% 2. Generate synthetic series, gaussian i.i.d. for example: 
testsig=300;
testmu=800;
nyr=100;
ytest=randn(nyr,1)*testsig + testmu; %(synthetic) 
drsti=round(0.4*nyr);dreni=round(0.55*nyr); %bad drought start
ytest(drsti:dreni)=0.5*ytest(drsti:dreni); %add a bad drought
ytest(ytest<0)=0; %check all data is positive
t=1901:1901+nyr-1; %years of ytest (synthetic)
varnm='Streamflow Volume (GL)'; %meaningful name and units (synthetic)
sitenm='Synthetic test site'; %site name

%% 3. Calculate and plot the TPM-derived ANEP curves 

%Set some parameters of TPM analysis:
LenSum=[3;5;10]; %m-year droughts to explore
NumClass=20; %number of TPM classes 
annualflag=true; %annual data
nomethod='RNOM'; %RNOM method for non-overlapping sum calculation
cbtype='DryEndFocus'; %produce distribution which focuses more on the dry end

%Set up the plot: 
f1=figure('position',[50 400 1300 900],'color','white');
ylm=[0 round(10*testmu,-3)];
probxlim=[0.01 0.6];
p_label = [0.01 0.05 0.1 0.25 0.5];
xtuse=norminv(p_label);
nwd=1; %number of worst obs droughts to show on plots
lw1=2;ms1=8; %plot settings

%Loop over LenSum values: 
for mm=1:length(LenSum)

  %Calculate the ANEP distributions:  
  [DroughtVol,AnnNEP,ClassBound,tpm,sumvolnos,sumvolnost,sumvolos,sumvolost]=...
    tpmdrought(ytest,LenSum(mm),NumClass,t,annualflag,nomethod,cbtype);
  
  % Plot the ANEP curve:   
  plAnnNEPi=AnnNEP<0.5; %only interested in the dry end
  pl(mm)=plot(norminv(AnnNEP(plAnnNEPi)),DroughtVol(plAnnNEPi),'o-',...
    'LineWidth',lw1,'markersize',ms1); hold all;
  leg{mm}=[num2str(LenSum(mm)) '-yr']; %legend entry
  set(gca,'xlim',norminv(probxlim)) %set x-axis up 
  xticks(xtuse); set(gca,'XTickLabelMode','Manual','XTickLabel',p_label,'XDir','Reverse')
  ylabel(varnm);xlabel('Annualised Non-Exceedence Probability (ANEP)');grid on;  
  yrtxt=['(' num2str(t(1)) '-' num2str(t(end)) ')'];
  title(['Annualised non-exceedence probabilities for multi-year sums for ' sitenm ' ' yrtxt])
  ylim(ylm); set(gca,'YTickLabelMode','Auto'); grid minor;
  legend(flip(pl),flip(leg))

end
 
%%  4. Plot the syntehtic time series for reference
figure('position',[1365  860 800 400]);
plot(t,ytest);title(sitenm);ylabel(varnm)



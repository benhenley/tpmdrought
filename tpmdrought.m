function [DroughtVol,CumNEP,AnnNEP,ClassBound,tpm,sumvolnos,sumvolnost,...
  sumvolos,sumvolost]=tpmdrought(vol,LenSum,NumClass,t,annualflag,...
  nomethod,cbtype)
% Code by Ben Henley: @benhenley, www.benhenley.net
% tpmdrought: Calculates exceedance probabilities of multi-timestep totals
%
% SYNTAX: 
% [DroughtVol,CumNEP,AnnNEP,ClassBound,tpm,sumvolnos,sumvolnost,...
%  sumvolos,sumvolost]=tpmdrought(vol,LenSum,NumClass,t,annualflag,...
%  nomethod,cbtype)
%
% INPUTS:  
% vol: input data timeseries as a [numm x 1] vector (e.g. streamflow, rainfall)
% LenSum: number of timesteps (e.g. years/months) for LenSum-yr or LenSum-month droughts 
% NumClass: number of states or classes used in the TPM calculation, 
%   equal to the number of cumulative probability values output
% numm: length of vol
% t: list of years or months of vol [numm x 1], e.g. 1930:1:1990 or 1:99 
% annualflag: true if annual data, false if monthly, (currently no other option).
% nomethod: char string which must match one of the following: 
%   'SNOS' option calculates sequential non-overlapping sums (SNOS) 
%     sequence of LenSum totals from the start to the end of the timeseries
%   'RNOM' option calculates sums using ranked non-overlapping minima (RNOM)
%     method, which finds the lowest non-overlapping sum in the series, then 
%     removes that sum from the time series and repeats until all LenSum sums are gone
% cbtype: the class bound setting, either 'DryEndFocus' or 'EqNumEvPerClass'%   
%   'EqNumEvPerClass' creates classes with an equal number of events per class
%   'DryEndFocus' creates classes which allow fewer events in the dry end of the distribution
%
% OUTPUTS: 
% DroughtVol: Volume to be plotted against AnnNEP, upper ClassBound (states) of TPM [NumClass x 1] 
% CumNEP: Cumulative Non-Exceedance probabilities of multi-year (timestep) droughts [NumClass x 1]
% AnnNEP: Annualised Non-Exceedance probabilities of multi-year (timestep) droughts [NumClass x 1]
% ClassBound: Upper boundaries of the classes of the TPM [NumClass x 1] 
% tpm: steady state transition probability matrix [NumClass x NumClass]
% sumvolno: non-overlapping totals using nomethod [numno x 1]
% sumvolnot: list of years or months of starting point of non-overlapping totals  
% sumvolos: overlapping totals sorted in ascending order
% sumvolost: list of years or months of starting point of overlapping totals  
%
% Code history: 
% v1.00: Initial commit, TPM code adapted by Ben Henley primarily from 
% Rory Nathan's fortran code; with Tom McMahon's generous assistance
%
% Disclaimer:
% This software is provided “as is”, without warranty of any kind, express 
% or implied, including but not limited to the warranties of merchantability, 
% fitness for a particular purpose and noninfringement. In no event shall 
% the authors or copyright holders be liable for any claim, damages or other 
% liability, whether in an action of contract, tort or otherwise, arising 
% from, out of or in connection with the software or the use or other dealings 
% in the software.

%1. Check inputs: ----
assert(isvector(vol),'tpmdrought:vol is not a vector');
assert(all(vol>=0),'tpmdrought:vol were not all zero')
numm=length(vol); assert(length(t)==numm,'tpmdrought:vol and t not same length')

%2. Initialisation: ----
 
  %Get overlapping and non-overlapping sums:  
  numod=numm-LenSum+1; 
  [sumvolnos,numno,sumvolnost]=SumNolap(LenSum,numm,numod,t,vol,nomethod);
  [sumvolo]=SumOlap(LenSum,numod,vol); 
  [sumvolos,indx]=sort(sumvolo);
  sumvolost=t(indx);
  sumvoluse=sumvolo;
  
  %Get Class bounds: (cbtype='DryEndFocus' or 'EqNumEvPerClass') 
  switch cbtype
    case 'EqNumEvPerClass' %equal number of events per class
    ClassBound=prctile(sumvoluse,100.*linspace(0,1,NumClass+1)');
    ClassBound(1)=0;
    ClassBound(end)=Inf;
    
    case 'DryEndFocus'      
      ClassBound=nan(NumClass+1,1);
      ClassBound(1)=0;
      NumEventLowestClass=max(5,length(sumvoluse)/100); %tinker with the number of events in the lowest class
      cbmin = NumEventLowestClass/length(sumvoluse); %one fifth up the distribution
      cbmax = 1.0-2.0/NumClass;
      ClassBound(2:NumClass)=prctile(sumvoluse,100.*linspace(cbmin,cbmax,NumClass-1)');
      ClassBound(end)=Inf;      
    otherwise %users can develop their own class bound methods 
      error('Invalid class bound type in tpmdrought (use DryEndFocus or EqNumEvPerClass)')      
  end  
    
  %Calculate TPM:
  tpm=calctpm(sumvoluse,ClassBound,NumClass);

  %Calculate volumes and probabilities quantiles
  CumNEP=nan(NumClass,1);AnnNEP=nan(NumClass,1);
  if annualflag
    adj=LenSum;
  else %monthly data
    adj=LenSum/12;
  end
   
  for ii=1:NumClass    
    CumNEP(ii)=sum(tpm(1:ii,1));
    AnnNEP(ii)=CumNEP(ii)./adj;
  end
  DroughtVol=ClassBound(2:end);  
  
end 

%% Function for Overlapping Sums
function [sumvol]=SumOlap(LenSum,numod,vol) 
  sumvol=nan(numod,1);
  for ii=1:numod
    sumvol(ii)=sum(vol(ii:ii+LenSum-1));
  end  
end

%% Function for Non-overlapping Sums
function [sumvols,numno,monrs]=SumNolap(LenSum,numm,numod,mon,vol,nomethod)

  sumvol=nan(numod,1);
  monr=nan(numod,1); 
  
  switch nomethod %switch on the method to extract non-overlapping sums 
    
    case 'SNOS'
  %The sequential non overlapping sums (SNOS) method extracts the 
  % non-overlapping sums end-to-end on each other as a sequence directly from the timeseries
    nn=1;ii=1; 
    while ii<=numod
      sumvol(nn)=sum(vol(ii:ii+LenSum-1));
      monr(nn)=mon(ii);
      ii=ii+LenSum;nn=nn+1;
    end
    numno=nn-1; sumvol=sumvol(1:numno,:); monr=monr(1:numno,:);
    [sumvols,inds]=sort(sumvol);
    monrs=monr(inds);    
    
    case 'RNOM'  
    % This method finds the lowest non-overlapping sum in the series, then 
    % removes that sum from the time series and repeats until all LenSum 
    % sums are gone 
    iflg=ones(numm,1); 
    ifound=1;
    n=1; 
    while ifound==1  
      smin=realmax;
      ifound=0;
      for ii=1:numod
        sumflg=sum(iflg(ii:ii+LenSum-1));
        if sumflg < LenSum 
          continue; 
        end
        ifound = 1;
        ss=sum(vol(ii:ii+LenSum-1));
        if ss < smin
          smin = ss;
          imin = ii;
        end
      end    
      if ifound==0
        break
      end
      sumvol(n)=smin;
      monr(n)=mon(imin);
      iflg(imin:imin+LenSum-1)=0;
      n=n+1;      
    end
    numno=n-1; sumvols=sumvol(1:numno,:); monrs=monr(1:numno,:);

    otherwise
      error(['nydrought: Invalid method for extracting ' ...
      'non-overlapping sums, choose nomethod as either ''SNOS'' or ''RNOM'''])
      
  end
  
end

%% Function for TPM Calc
function tpm=calctpm(sumvol,ClassBound,NumClass)

  tpm=nan(NumClass,NumClass);
  nc=zeros(NumClass,NumClass);      
  sumvol=makerow(sumvol); %strictly a row vector
  ClassBound=makecol(ClassBound); %strictly a col vector
  cl=sum(sumvol>=ClassBound); %get the list of classes using a simple linear algebra operation
  
  for ii=1:length(cl)-1 %use the class list to index the class counts matrix   
    nc(cl(ii),cl(ii+1))=nc(cl(ii),cl(ii+1))+1;
  end      

  for ii=1:NumClass
    xn=sum(nc(ii,1:NumClass));
    for jj=1:NumClass
      tpm(ii,jj)=nc(ii,jj)/xn;
    end
  end
  tpm=transpose(tpm);
  tpm=ForceUnity(tpm,NumClass);

  %calculate steady state TPM
  for ii=1:20  
    tpm=tpm^2;
    tpm=ForceUnity(tpm,NumClass);
  end

end

%% Function to Force unity matrix
function tpm=ForceUnity(tpm,NumClass)
  for jj=1:NumClass
    spx=sum(tpm(1:NumClass,jj));
    for ii=1:NumClass
      tpm(ii,jj)=tpm(ii,jj)/spx;
    end
  end 
end    

%END ----    

function [PDSI,X1,X2,X3,Pe,montho]=PDSIr(Z);
nmonths=0;
X1=zeros(size(Z));X2=X1;X3=X1;Uw=X1;Ud=X1;Pe=X1;PDSI=X1;Ze=X1;
% start off not in dry or wet spell

XX=Z(1)/3;
if XX>0 X1(1)=XX;else X2(1)=XX;end
if abs(XX)>=1 X3(1)=XX;nmonths=1;end
PDSI(1)=XX;


% loop over time

for i=2:length(Z)
XX=Z(i)/3;
Ud(i)=XX*3-.15;
Uw(i)=XX*3+.15;
if XX>0 
 if X2(i-1)<0 X2(i)=X2(i-1)*.897+XX;else X2(i)=XX;end; if X2(i-1)<-1 & X3(i-1)<-1 X2(i)=0;end
 if X1(i-1)>0 X1(i)=X1(i-1)*.897+XX;else X1(i)=XX;end; if X2(i-1)<-1 & X3(i-1)<-1 X2(i)=0;end
  if (X3(i-1))~=0
   if X3(i-1)<0
    Ze(i)=-2.691*X3(i-1)-1.5;
    Pe(i)=gettingout(Uw,Ze,nmonths+1,i,1);
    if Pe(i)==100 nmonths=0; X3(i)=0; if X1(i)>1 X3(i)=X1(i);nmonths=1; end;end
    if Pe(i)==0 X1(i)=0;X2(i)=0;end
   else
    Ze(i)=-2.691*X3(i-1)+1.5;
    Pe(i)=gettingout(Ud,Ze,nmonths+1,i,0);
    if Pe(i)==100 nmonths=0; X3(i)=0; if X1(i)>1 X3(i)=X1(i);nmonths=1; end;end
    if Pe(i)==0 X1(i)=0;X2(i)=0;end    
   end
  else
   if X1(i)>0.5 X3(i)=X1(i);nmonths=1;end
  end
else
 if X2(i-1)<0 X2(i)=X2(i-1)*.897+XX;else X2(i)=XX;end
 if X1(i-1)>0 X1(i)=X1(i-1)*.897+XX;else X1(i)=XX;end
 if (X3(i-1))~=0
  if X3(i-1)>0
   Ze(i)=-2.691*X3(i-1)+1.5;
   Pe(i)=gettingout(Ud,Ze,nmonths+1,i,0);
   if Pe(i)==100 nmonths=0; X3(i)=0;if X2(i)<-1 X3(i)=X2(i); nmonths=1;end;end
   if Pe(i)==0  X1(i)=0;X2(i)=0;end
   else
    Ze(i)=-2.691*X3(i-1)-1.5;
    Pe(i)=gettingout(Uw,Ze,nmonths+1,i,1);
    if Pe(i)==100 nmonths=0; X3(i)=0; if X2(i)>1 X3(i)=X2(i);nmonths=1; end;end
    if Pe(i)==0 X1(i)=0;X2(i)=0;end    
  end
 else
  if X2(i)<-0.5 X3(i)=X2(i);nmonths=1;end
 end
end
if Pe(i)<100
 if (X3(i-1))~=0 X3(i)=X3(i-1)*.897+XX;nmonths=nmonths+1;end
end
% decide what PDSI is
montho(i)=nmonths;
% option 1 no drought going on or being established, nmonths=0
%PDSI(i)=PDSI(i)*.897+XX;
if X1(i)<0 X1(i)=0;end 
if X2(i)>0 X2(i)=0;end

if nmonths==0
 if X3(i)==0 
   if X1(i)>-X2(i) PDSI(i)=X1(i);else PDSI(i)=X2(i);end
   PDSI(i)=XX;
 end
else
% if in wet/dry spell and it has not ended PDSI=X3
%  if abs(X3(i))<1
%      if X1(i)>-X2(i) PDSI(i)=X1(i);else PDSI(i)=X2(i);end
%  else
% if nmonths>1
  PDSI(i)=X3(i);
% end

 if PDSI(i) >=1 & nmonths>1 X1(i)=0; end
 if PDSI(i) <=-1 & nmonths>1 X2(i)=0; end
end
end% if beginning wet/dry spell PDSI=X1 or X2

  


% if beginning wet/dry spell PDSI=X1 or X2



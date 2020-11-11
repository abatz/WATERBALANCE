function [AET,DEF,RUNOFF,SNOWS,SOILS,RUNOFFSNOW]=simplehydromodel(TMEAN,PPT,PET,AWC,soil,snow);
% PPT and PET should be in mm
% data can be daily/monthly
% simple 1-bucket soil-runoff model
DR=NaN;

AET=single(NaN*ones(size(PPT)));DEF=AET;RUNOFF=AET;IN=AET;RUNOFFSNOW=AET;
% start with empty soil bucket
snowstorage=snow;
for j=1:size(TMEAN,2)
TMN=TMEAN(:,j);
% initialize drainsnow and snowdrink
drainsoil=0*TMN;
snowdrink=0*TMN;

MF=1-runsnow(TMN+273.15,1);
SNOW=(1-MF).*PPT(:,j);RAIN=MF.*PPT(:,j);
MELT=MF.*(SNOW+snowstorage);
INPUT=RAIN+MELT;
FRACTRAIN=max(0,RAIN./INPUT);
f=find(isinf(FRACTRAIN==1));FRACTRAIN(f)=1;
f=find(isnan(FRACTRAIN==1));FRACTRAIN(f)=1;
FRACTRAIN=min(1,FRACTRAIN);
% assume that 5% of input goes to direct runoff
extrarun=INPUT*0.05;
INPUT=INPUT*.95;

option=1;

% how much of this is from snow
% option 1, prioritize all from snow

if option==1
RUNOFFSNOW(:,j)=min(extrarun,MELT);
RRAIN=extrarun-RUNOFFSNOW(:,j);
%MELT=max(MELT-extrarun,0);
RAININPUT=RAIN-RRAIN;
end
    
    
snowstorage=(1-MF).*(SNOW+snowstorage);
deltasoil=INPUT-PET(:,j);
if option==1 excessafterliquid=max(0,RAININPUT-PET(:,j));end
f1=find(deltasoil<0); % if water needs to be taken from soil
f2=find(snowstorage(f1)>0 & snowstorage(f1)>-deltasoil(f1)); %if all can be taken from snow layer
  snowstorage(f1(f2))=snowstorage(f1(f2))+deltasoil(f1(f2));
  snowdrink(f1(f2))=-deltasoil(f1(f2));
  deltasoil(f1(f2))=0;
f3=find(snowstorage(f1)>0 & snowstorage(f1)<-deltasoil(f1)); % if soil moisture needs to be drained as well
  deltasoil(f1(f3))=deltasoil(f1(f3))+snowstorage(f1(f3));
  snowdrink(f1(f3))=snowstorage(f1(f3));
  snowstorage(f1(f3))=0;

% after eating snowpack, what else is needed from soil moisture
f1=find(deltasoil<0);f2=find(deltasoil>0);
    ff=find(-deltasoil>soil);
% can not drain more than is in soil!!!
    deltasoil(ff)=-soil(ff);
    drainsoil(f1)=deltasoil(f1).*(1-exp(-soil(f1)./AWC(f1)));
    drainsoil(f2)=0;

demand=PET(:,j);
supply=INPUT+snowdrink; % supply without mining soil water


f=find(demand>=supply);
f1=find(demand<supply);
 AET(f,j)=supply(f)-drainsoil(f);
 DEF(f,j)=PET(f,j)-AET(f,j);
 RUNOFF(f,j)=0;
 soil(f)=soil(f)+drainsoil(f);

 AET(f1,j)=PET(f1,j);
 DEF(f1,j)=0;
f2=find(soil(f1)+deltasoil(f1)>AWC(f1));
f3=find(soil(f1)+deltasoil(f1)<=AWC(f1));
if option==1
  excess=max(0,soil+deltasoil-AWC);
  excessrainonly=max(0,soil+excessafterliquid-AWC);
  RUNOFF(f1(f2),j)=excess(f1(f2));
  RUNOFFSNOW(f1(f2),j)=RUNOFFSNOW(f1(f2),j)+(excess(f1(f2))-excessrainonly(f1(f2)));
else
  excess=max(0,soil+deltasoil-AWC);
  RUNOFF(f1(f2),j)=excess(f1(f2));
  RUNOFFSNOW(:,j)=excess.*(1-FRACTRAIN);
  RUNOFFSNOW(:,j)=RUNOFFSNOW(:,j)+(1-FRACTRAIN).*extrarun;
end
  soil(f1(f2))=AWC(f1(f2));

   soil(f1(f3))=soil(f1(f3))+deltasoil(f1(f3));
   RUNOFF(f1(f3),j)=0;
SOILS(:,j)=soil;
SNOWS(:,j)=snowstorage;
RUNOFF(:,j)=RUNOFF(:,j)+extrarun;
end;

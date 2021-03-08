function [SNOWS]=snowrun(TMEAN,PPT,snow);
% PPT and PET should be in mm
% data can be daily/monthly
% simple 1-bucket soil-runoff model
DR=NaN;

snowstorage=snow;
for j=1:size(TMEAN,2)
TMN=TMEAN(:,j);
% initialize drainsnow and snowdrink
drainsoil=0*TMN;
snowdrink=0*TMN;

MF=1-runsnow(TMN+273.15,1);
SNOW=(1-MF).*PPT(:,j);RAIN=MF.*PPT(:,j);
MELT=MF.*(SNOW+snowstorage);
SNOWS(:,j)=(1-MF).*(SNOW+snowstorage);
end;

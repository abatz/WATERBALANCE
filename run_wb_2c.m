
% if you want to apply CO2 resistance correction for elevated CO2
% Kruijt, B., Witte, J.P.M., Jacobs, C.M. and Kroon, T., 2008. Effects of rising atmospheric CO2 on evapotranspiration and soil moisture: A practical approach for the Netherlands. Journal of Hydrology, 349(3-4), pp.257-267.
if isvar(degrees)
if degrees==2
 petdata=petdata*.978;
 dd='2';
else
 petdata=petdata*.948;
 dd='4';
end
end

% load soil water holding capacity layer - 
soilt=single(soilt);

% loop over years to load in tmax/tmin/pet/ppt data
for yr=1:31
 tmean=tmaxdata/2+tmindata/2;
 clear tmaxdata tmindata dataout
 tmean=shiftdim(tmean,2);
 pptdata=shiftdim(pptdata,2);
 petdata=shiftdim(petdata,2);

% f refers to non-empty cells if needed

 tmean=tmean(:,f);
 pptdata=pptdata(:,f);
 petdata=single(petdata(:,f));
 soilt=soilt(f);

% carry forward soil and snow from previous december or spinup if needed; otherwise start with some default soil water and snow

 [AET.data,DEF.data,RO.data,SNOW.data,SOIL.data,ROSNOW.data]=hydro_tax_ro(tmean',pptdata',petdata',soilt',soillast',snowlast');
soillast=soildata(12,:);
snowlast=snowdata(12,:);

aetdata=round(aetdata,1);
defdata=round(defdata,1);
runoffdata=round(runoffdata,1);
snowdata=round(snowdata,1);
soildata=round(soildata,1);
rosnowdata=round(rosnowdata,1);

save([dirr,'wb_',dd,'cyue_',num2str(1984+yr)],'-v7.3','aetdata','defdata','runoffdata','snowdata','soildata','rosnowdata');
clear *data*
end

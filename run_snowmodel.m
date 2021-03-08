% loop over years to load in tmax/tmin/pet/ppt data
for yr=1:31
 tmean=tmaxdata/2+tmindata/2;
 clear tmaxdata tmindata dataout
 tmean=shiftdim(tmean,2);
 pptdata=shiftdim(pptdata,2);

% f refers to non-empty cells if needed

 tmean=tmean(:,f);
 pptdata=pptdata(:,f);

% carry forward snow from previous december or spinup if needed; otherwise start with some default snow

 [ROSNOW.data]=snowrun(tmean',pptdata',snowlast');
snowlast=snowdata(12,:);

snowdata=round(snowdata,1);


save([dirr,'swe_',num2str(1984+yr)],'-v7.3','snowdata');
clear *data*
end

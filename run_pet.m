function [petdata]=run_pet(vap,tmax,tmin,srad,ws,el);

% approximate surface pressure from elevation
P=1013.25*((293-0.0065*el)/293).^5.26;

 vpn=calcVP(tmin+273.15,repmat(P,[1 1 12]));
 vpx=calcVP(tmax+273.15,repmat(P,[1 1 12]));
 vpn=vpn/1000;vpx=vpx/1000;
 vpd=vpn/2+vpx/2-vap;
 vpd(vpd<0)=0;
 clear vap vpx vpn clear vap
 tmax=shiftdim(tmax,2);
 tmin=shiftdim(tmin,2);
 ws=shiftdim(ws,2);
 srad=shiftdim(srad,2);
 vpd=shiftdim(vpd,2);
 petdata=NaN*ones(12,4320,8640);
 % exclude ocean points
 f=find(~isnan(tmax(1,:))==1);
 vpd=vpd(:,f);
 tmax=tmax(:,f);
 tmin=tmin(:,f);
 srad=srad(:,f);
 ws=ws(:,f);
 [ET]=monthlyPETvpd(srad'/86.4, tmax',tmin', ws',lat(f)',el(f)',0.23,vpd');
 petdata(:,f)=ET';
 
 % finally we restrict PET that occurs @ low temperatures for water balance
 % purposes; this may or not be needed depending on application
 t=1-runsnow(tmax/2+tmin/2,ones(size(tmax)));
 petdata(:,f)=ET'.*t';
 

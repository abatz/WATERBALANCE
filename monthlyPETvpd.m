function [ET]=monthlyPET(radiation, tmax,tmin, wind,Lat,Z,albedo,vpd);
% Calculates potential evapotranspiration using Pennman-Montieth equation
% source: Allen et al. 1998
% Input: (radiation in MJ/m2/d), tmax,tmin in C, wind in m/s, month
% lasttmean is temperature difference from last month
% Lon, Lat, Z 
% Uses function nanmean

% convert radiation from W/m2 to MJ/d 

nn=size(radiation);
if ndims(radiation)==3
 s1=size(radiation);
 radiation=reshape(shiftdim(radiation,2),12,s1(1)*s1(2))';
 tmax=reshape(shiftdim(tmax,2),12,s1(1)*s1(2))';
 tmin=reshape(shiftdim(tmin,2),12,s1(1)*s1(2))';
 wind=reshape(shiftdim(wind,2),12,s1(1)*s1(2))';
 vpd=reshape(shiftdim(vpd,2),12,s1(1)*s1(2))';
 Lat=Lat(:);Z=Z(:);
end

radiation=radiation*.0864;

daysinmonth=[31    28    31    30    31    30    31    31    30    31    30    31];
d2=[31 59 90 120 151 181 212 243 273 304 334 365];
d1=[1 32 60 91 121 152 182 213 244 274 305 335];

% calculate change in temperature from month to month
tmean=(tmax+tmin)/2;

% calculate change in temperature from month to month
lasttmean=NaN*ones(size(tmax));
tmean=(tmax+tmin)/2;
for i=2:size(tmax,2)
    lasttmean(:,i)=tmean(:,i)-tmean(:,i-1);
end
lasttmean(:,1)=tmean(:,1)-tmean(:,size(tmean,2));

% wind adjustment to 2m from 10m output, can be changed for ET at
% elevations other than 2m, but need to scale using logarithmic increase in
% wind speed with elevation

wind=wind*(4.87/log(67*10-5.42));

% Saturation vapor pressure , assume sat vap pressure is TMIN-2
es1 = 0.6108 * exp(tmin * 17.27./ (tmin + 237.3));
es2 = 0.6108 * exp(tmax * 17.27./ (tmax + 237.3));
es = es1/2.+es2/2.;
%ea = 0.6108 * exp(tdew * 17.27./ (tdew + 237.3));
%vpd=es-ea;
ea=es-vpd;
g=find(vpd<0);vpd(g)=0;
% VPD - Vapor pressure deficit 
%ea=es-vpd;
%VPD = es - ea; % (kPa) 
% DEL - Slope of the saturation vapor pressure vs. air temperature curve at the average hourly air temperature 
DEL = (4098 * es)./(tmean + 237.3).^2;


% Barometric pressure 
P = 101.3*((293-0.0065*Z)/293).^5.26; % CIMIS
lambda=2.501-2.361e-3*tmean;
% GAM - Psychrometer constant (kPa C-1) \

GAM = 0.00163*repmat(P,[1 12])./lambda;
%GAM2=6.65e-4.*repmat(P',[1 12]);
%GAM = 0.000646 * (1 + 0.000946*tmean) .* repmat(P',[1 12]); % CIMIS
% W - Weighting function 
%W = DEL./(DEL + GAM);



GSC = 0.082; % MJ m -2 min-1 (solar constant)

phi = pi*Lat/180;

for doy=1:12
% should average this over a month
% Calculate potential max solar radiation or clear sky radiation
% assumed to be 75% TOA shortwave radiation or cloudless day, FAO, 1998
clear Rso
for i=1:daysinmonth(doy)
    DoY=d1(doy)-1+i;
    dr = 1+0.033*cos(2*pi/365 * DoY); 
    delta = 0.409 * sin(2*pi/365*DoY-1.39);
    omegas = acos(-tan(phi).*tan(delta));
    Ra = 24*60/pi.*GSC.*dr .* ( omegas .*sin(phi).*sin(delta) +cos(phi).*cos(delta).*sin(omegas) ); % FAO daily
    Rso(:,i) = Ra .* (0.75+2e-5*Z);
end
% incomming solar radiation has already been corrected for macroscale
% albedo, may need to calibrate for actual albedo

Rso=real(mean(Rso,2));
% incomming solar radiation has already been corrected for macroscale
% albedo, may need to calibrate for actual albedo

f=find(Rso<0 | isnan(Rso)==1);
Rso(f)=0;


% lr out, Rso>=radiation
% radfraction is a measure of relative shortwave radiation, or % of
% possible radiation, needs to be less than 1
radfract=radiation(:,doy)./Rso;
f=find(radfract>1);radfract(f)=1;
f=find(isinf(radfract)==1);radfract(f)=1;
f=find(isnan(radfract)==1);radfract(f)=1;
longw=4.903e-9*((tmax(:,doy)+273.15).^4+(tmin(:,doy)+273.15).^4)/2.*(.34-.14*sqrt(ea(:,doy))).*(1.35*radfract-.35);
%[nmean(longw(:))]
% netrad
netrad=real(radiation(:,doy).*(1-albedo)-longw)*daysinmonth(doy);
%nmean(netrad(:))
%soil heat flux
G=0.14*(lasttmean(:,doy))*daysinmonth(doy);

%ET
TERM1_NUMERATOR = .408*DEL(:,doy).*(netrad-G);
TERM2_NUMERATOR=GAM(:,doy).*wind(:,doy).*vpd(:,doy).*900./(273.15+tmean(:,doy))*daysinmonth(doy);
DENOMINATOR=DEL(:,doy)+GAM(:,doy).*(1+0.34.*wind(:,doy));
TERM1=TERM1_NUMERATOR./DENOMINATOR;
TERM2=TERM2_NUMERATOR./DENOMINATOR;
ET(:,doy) = TERM1+TERM2;
%[nmean(TERM1) nmean(TERM2)]
end
f=find(ET<0);ET(f)=0;
ET=reshape(ET,nn);

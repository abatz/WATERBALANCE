function [vpd]=calcVPD(Ta,P);
% SATVAP: computes saturation vapor pressure
% q=satvap(Ta) computes the vapor pressure at satuation at air
% temperature Ta (deg C). From Gill (1982), Atmos-Ocean Dynamics, 606.
%
%    INPUT:   Ta- air temperature  [C]
%             p - pressure (optional)  [mb]
%
%    OUTPUT:  q  - saturation vapour pressure  [mb]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
% 8/27/98: version 1.1 (corrected by RP)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
P=1013.25;
end
Ta=Ta-273.16;

ew=power(10,((0.7859+0.03477*Ta)./(1+0.00412*Ta)));

fw=1 + 1e-6*P.*(4.5+0.0006.*Ta.^2);
vpd=fw.*ew*100;



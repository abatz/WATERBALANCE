function [snow]=runsnow(temp,precip);
a=-48.2292;
b=0.7205;
c=1.1662;
d=1.0223;
temp=temp-273.15;
snowf=a*(tanh(b*(temp-c))-d);
f=find(temp<-2);snowf(f)=100;
f=find(temp>6.5);snowf(f)=0;
snow=single(snowf./100.*precip);


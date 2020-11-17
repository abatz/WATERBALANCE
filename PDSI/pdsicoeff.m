function[alp,bet,gam,del]=ja_cc_0801(PESUM,ETSUM,RSUM,PRSUM,SPSUM,ROSUM,PLSUM,TLSUM)
f1=find(PESUM~=0);
if ~isempty(f1)
alp(f1)=ETSUM(f1)./PESUM(f1);
end
f2=setxor(1:length(PESUM),f1);
if ~isempty(f2)
f3=find(ETSUM(f2)==0);
if ~isempty(f3)
alp(f2(f3))=1.0;
end
f4=find(ETSUM(f2)~=0);
if ~isempty(f4)
alp(f2(f4))=0.0;
end;
end
%Beta Calcuation
f5=find(PRSUM~=0);
if ~isempty(f5)
bet(f5)=RSUM(f5)./PRSUM(f5);
end
f6=setxor(1:length(PRSUM),f5);
if ~isempty(f6)
f7=find(RSUM(f6)==0);
if ~isempty(f7)
bet(f6(f7))=1.0;
end
f8=find(RSUM(f6)~=0);
if ~isempty(f8)
bet(f6(f8))=0.0;
end;
end
%Gamma Calcuation
f9=find(SPSUM~=0);
if ~isempty(f9)
gam(f9)=ROSUM(f9)./SPSUM(f9);
end


f10=setxor(1:length(SPSUM),f9);
if ~isempty(f10)
f11=find(ROSUM(f10)==0);

if ~isempty(f11)
gam(f10(f11))=1.0;
end
f12=find(ROSUM(f10)~=0);
if ~isempty(f12)
gam(f10(f12))=0.0;
end;
end
%Delta Calcuation
f13=find(PLSUM~=0);
if ~isempty(f13)
del(f13)=TLSUM(f13)./PLSUM(f13);
end
f14=setxor(1:length(PLSUM),f13);
if ~isempty(f14)
del(f14)=0.0;
end;

%alp=ETSUM./PESUM;
%bet=RSUM./PRSUM;
%del=TLSUM./PLSUM;
%gam=ROSUM./SPSUM;


alp=alp';
bet=bet';
gam=gam';
del=del';
return


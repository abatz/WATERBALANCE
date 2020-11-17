function[PDSI,SOIL,Z_1]=calc_pdsi(PET,ppt1,WCTOP,WCBOT,ss2)
%INITIALIZE VARS
WCTOT=WCBOT+WCTOP;
SS=WCTOP;
SU=WCBOT;
SPSUM=0;SPSUM=repmat(SPSUM,[size(WCTOT,1) 12]);
PLSUM=0;PLSUM=repmat(PLSUM,[size(WCTOT,1) 12]);
PRSUM=0;PRSUM=repmat(PRSUM,[size(WCTOT,1) 12]);
RSUM=0;RSUM=repmat(RSUM,[size(WCTOT,1) 12]);
TLSUM=0;TLSUM=repmat(TLSUM,[size(WCTOT,1) 12]);
ETSUM=0;ETSUM=repmat(ETSUM,[size(WCTOT,1) 12]);
ROSUM=0;ROSUM=repmat(ROSUM,[size(WCTOT,1) 12]);
DBAR=0;DBAR=repmat(DBAR,[size(WCTOT,1) 1]);
SABSD=0;SABSD=repmat(SABSD,[size(WCTOT,1) 12]);
AK=0;AK=repmat(AK,[size(WCTOT,1) 12]);
AKHAT=0;AKHAT=repmat(AKHAT,[size(WCTOT,1) 12]);
SWTD=0;SWTD=repmat(SWTD,[size(WCTOT,1) 12]);
SP=SS+SU;
PSUM=nansum(ppt1,3);
PESUM=nansum(PET,3);

ssz=size(ppt1,3);
if nargin==4
ss2=ssz;
end
% fake ten year spinup
PET(:,:,11:ssz+10)=PET;PET(:,:,1:10)=repmat(nanmean(PET,3),[1 1 10]);
ppt1(:,:,11:ssz+10)=ppt1;ppt1(:,:,1:10)=repmat(nanmean(ppt1,3),[1 1 10]);

%BEGIN REAL STUFF
for yr=1:(size(ppt1,3))%1:200
%if ~isnan(PET(1,1,yr))
    for mo=1:12
        PR=WCTOT-SP;
        [PL,ET,TL,RO,R,SSS,SSU]=pdsi_wb(PET(:,mo,yr),ppt1(:,mo,yr),WCBOT(:),WCTOP(:),WCTOT(:),SS(:),SU(:),SP(:));
        SS=SSS;
        SU=SSU;
        SP=SS+SU;
        spdat(:,mo,yr)=single(SP);
        pldat(:,mo,yr)=single(PL);
        prdat(:,mo,yr)=single(PR);
        rdat(:,mo,yr)=single(R);
        tldat(:,mo,yr)=single(TL);
        etdat(:,mo,yr)=single(ET);
        rodat(:,mo,yr)=single(RO);
        sssdat(:,mo,yr)=single(SSS);
        ssudat(:,mo,yr)=single(SSU);
    end;
end;%clear SP PL PR R TL ET RO SS SU SSS SSU WCBOT WCTOP WCTOT yr mo
%end
c1=find(rdat>=prdat);rdat(c1)=prdat(c1);clear c1
c1=find(tldat>=pldat);tldat(c1)=pldat(c1);clear c1
ini=10+ss2(1):10+ss2(length(ss2));
SPSUM=nanmean(spdat(:,:,ini),3);
PLSUM=nanmean(pldat(:,:,ini),3);
PRSUM=nanmean(prdat(:,:,ini),3);
RSUM=nanmean(rdat(:,:,ini),3);
TLSUM=nanmean(tldat(:,:,ini),3);
ETSUM=nanmean(etdat(:,:,ini),3);
PESUM=nanmean(PET(:,:,ini),3);
ROSUM=nanmean(rodat(:,:,ini),3);
PSUM=nanmean(ppt1(:,:,ini),3);

%CAFEC
for mo=1:12
    [alp(:,mo),bet(:,mo),gam(:,mo),del(:,mo)]=pdsicoeff(PESUM(:,mo),ETSUM(:,mo),...
        RSUM(:,mo),PRSUM(:,mo),SPSUM(:,mo),ROSUM(:,mo),...
        PLSUM(:,mo),TLSUM(:,mo));
end;clear mo


TRAT=(PESUM+RSUM+ROSUM)./(PSUM+TLSUM);clear PESUM RSUM ROSUM PSUM TLSUM ETSUM PRSUM SPSUM PLSUM
PHAT=PET.*repmat(alp,[1 1 size(PET,3)])+prdat.*repmat(bet,[1 1 size(PET,3)])+spdat.*repmat(gam,[1 1 size(PET,3)])-pldat.*repmat(del,[1 1 size(PET,3)]);
DD=ppt1-PHAT;
SABSD=abs(DD);

DBAR=nanmean(SABSD,3);
AKHAT=1.5*log10((TRAT+2.8*25.4)./DBAR)+0.5;
f= AKHAT<0;AKHAT(f)=0;clear f
SWTD=sum(DBAR.*AKHAT,2);          

AK=17.67*25.4*AKHAT./repmat(SWTD,[1 12]);
Z_1=DD/25.4.*repmat(AK,[1 1 size(DD,3)]);
SOIL=reshape(sssdat+ssudat,size(sssdat,1),size(sssdat,2)*size(sssdat,3));

clear AK AKHAT DBAR DD PET SABSD SWTD alp bet del gam pldat ppt1 prdat rdat rodat sssdat ssudat tldat %etdat %spdat

% limit Z scores to +/- 16
f10= Z_1>16;Z_1(f10)=16;clear f10
f10= Z_1<-16;Z_1(f10)=-16;clear f10
PDSI=NaN*ones(size(Z_1,1),size(Z_1,2)*size(Z_1,3));
for i=1:size(PDSI,1)%parfor
    if max(Z_1(i,:))>0
        [PDSI(i,:),X1,X2,X3,Pe,montho]=PDSIr((Z_1(i,:)));
    end;
end;clear i
%matlabpool close
PDSI=PDSI(:,121:size(PDSI,2));
SOIL=SOIL(:,121:size(SOIL,2));
size(Z_1)
Z_1=Z_1(:,121:size(Z_1,2)*size(Z_1,3));
return

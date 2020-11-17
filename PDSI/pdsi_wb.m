function[PL,ET,TL,RO,R,SSS,SSU]=pdsi_wb(PET,P,WCBOT,WCTOP,WCTOT,SS,SU,SP)
%Created 06/18/11 by Jacob Wolf at the University of Idaho
%Recharge, Runoff, Residual Moisture, Loss to both SFC and Under Layers,
%depending on starting moisture content and vals of PPT and PET.
%from Palmer 1965.

PET=PET';
P=P';
WCBOT=WCBOT';
WCTOP=WCTOP';
WCTOT=WCTOT';
SS=SS';
SU=SU';
SP=SP';

% check on inverting matrices
% bottomfactor, loss

PR=WCTOT-SP;
PRS=WCTOP-SS;
PRU=WCBOT-SU;
PL=NaN*ones(size(P));ET=PL;TL=PL;RO=PL;R=PL;SSS=PL;SSU=PL;
f1=find(SS>=PET);
if ~isempty(f1)
    PL(f1)=PET(f1);
end
f2=find(SS<PET);
if ~isempty(f2)   
    straw=SU(f2)./WCTOT(f2);
    demand=PET(f2)-SS(f2);
    f5=find(demand>SU(f2));
    demand(f5)=SU(f2(f5));
    PL(f2)=demand.*straw+SS(f2);
    f4=find(PL(f2)>SP(f2));PL(f2(f4))=SP(f2(f4));clear straw demand
end
test2=find(P>=PET);
if ~isempty(test2)
% PPT exceeds PET
    ET(test2)=PET(test2);
    TL(test2)=0.0;
    test3=find((P(test2)-PET(test2))>(PRS(test2)));
    if ~isempty(test3)
    % PPT recharges under and upper layers
        RS(test2(test3))=PRS(test2(test3));
        SSS(test2(test3))=WCTOP(test2(test3));
        test4=find((P(test2(test3))-PET(test2(test3))-RS(test2(test3)))<(PRU(test2(test3))));
%        test4=find((P(test2(test3))-PET(test2(test3))-RS(test2(test3)))<(WCBOT(test2(test3))-SU(test2(test3))));
        if ~isempty(test4)
        % Both layers can take entire excess
            RU(test2(test3(test4)))=P(test2(test3(test4)))-PET(test2(test3(test4)))-RS(test2(test3(test4)));
            RO(test2(test3(test4)))=0.0;
        end
        test5=setxor(1:length(test2(test3)),test4);% Some runoff occurs
        if ~isempty(test5)
            RU(test2(test3(test5)))=WCBOT(test2(test3(test5)))-SU(test2(test3(test5)));
            RO(test2(test3(test5)))=P(test2(test3(test5)))-PET(test2(test3(test5)))-RS(test2(test3(test5)))-RU(test2(test3(test5)));
        end
        SSU(test2(test3))=SU(test2(test3))+RU(test2(test3));
        R(test2(test3))=RS(test2(test3))+RU(test2(test3));
    end    
    test6=setxor(1:length(test2),test3);% Only top layer is recharged
    if ~isempty(test6)
        R(test2(test6))=P(test2(test6))-PET(test2(test6));
        SSS(test2(test6))=SS(test2(test6))+P(test2(test6))-PET(test2(test6));
        SSU(test2(test6))=SU(test2(test6));
        RO(test2(test6))=0.0;
    end
end
f12=find(R>PR);R(f12)=PR(f12);

testa=setxor(1:length(P),test2);% Evaporation exceeds PPT
if ~isempty(testa)
    R(testa)=0.0;
    testb=find(SS(testa)>=(PET(testa)-P(testa)));
    if ~isempty(testb)
    % Evap from surface layer only
        SL(testa(testb))=PET(testa(testb))-P(testa(testb));
        SSS(testa(testb))=SS(testa(testb))-SL(testa(testb));
        UL(testa(testb))=0.0;
        SSU(testa(testb))=SU(testa(testb));
    end   
    testc=setxor(1:length(testa),testb);%Evap from both layers
    if ~isempty(testc)
        SL(testa(testc))=SS(testa(testc));
        SSS(testa(testc))=0.0;
    demand=PET(testa(testc))-P(testa(testc))-SL(testa(testc));

    f5=find(demand>SU(testa(testc)));
    demand(f5)=SU(testa(testc(f5)));
    straw=SU(testa(testc))./WCTOT(testa(testc));
    drainsoil=-demand.*straw;
    SSU(testa(testc))=SSU(testa(testc))+drainsoil;
    UL(testa(testc))=-drainsoil;
        f4=find(UL(testa(testc))>SU(testa(testc)));
        if ~isempty(f4)
            UL(testa(testc(f4)))=SU(testa(testc(f4)));
        end
        SSU(testa(testc))=SU(testa(testc))-UL(testa(testc));
    end  
    TL(testa)=SL(testa)+UL(testa);
    RO(testa)=0.0;
    ET(testa)=P(testa)+SL(testa)+UL(testa);
    f=find(PET(testa)<ET(testa));ET(testa(f))=PET(testa(f));
end

%ET(f)=PET(f);
ET=ET';
R=R';
RO=RO';
SSS=SSS';
SSU=SSU';
TL=TL';
PL=PL';
f=find(isnan(PET));
ET(f)=NaN;
R(f)=NaN;
RO(f)=NaN;
SSS(f)=NaN;
SSU(f)=NaN;
TL(f)=NaN;
PL(f)=NaN;
%f=find(ET>PET);ET(f)=PET(f);
%f=find(TL>PL);TL(f)=PL(f);
return

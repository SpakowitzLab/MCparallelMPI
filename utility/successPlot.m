path=sprintf('../ReptAndSwap7_19_24/');

repMIN=1;
repSkip=1;
repMAX=47;

nrep=length([repMIN:repSkip:repMAX]);

Chi=zeros(nrep,1);
CrankShaftSuc=zeros(nrep,1);
CrankShaftAmp=zeros(nrep,1);
SlideSuc=zeros(nrep,1);
PivotSuc=zeros(nrep,1);
FCRotSuc=zeros(nrep,1);
FCRotAmp=zeros(nrep,1);
FCTranSuc=zeros(nrep,1);
FCTranAmp=zeros(nrep,1);
SwapSuc=zeros(nrep,1);
ReptSuc=zeros(nrep,1);


for rep=[repMIN:repSkip:repMAX]
    fprintf('on %d of %d\n',rep,nrep)
    out=dlmread(strcat(path,sprintf('/out3v%d',rep)),'',1,0);
    out=out';
    CrankShaftSuc(rep)=out(5,end);
    CrankShaftAmp(rep)=out(4,end);
    SlideSuc(rep)=out(8,end);
    PivotSuc(rep)=out(11,end);
    FCRotSuc(rep)=out(17,end);
    FCRotAmp(rep)=out(16,end);
    FCTranSuc(rep)=out(20,end);
    FCTranAmp(rep)=out(19,end);
    SwapSuc(rep)=out(26,end);
    ReptSuc(rep)=out(28,end);
    out=dlmread(strcat(path,sprintf('/out1v%d',rep)),'',1,0);
    Chi(rep)=out(end,13);
end

figure(1)
hold on
plot(Chi,CrankShaftSuc,Chi,SlideSuc,Chi,PivotSuc,Chi,FCRotSuc,...
     Chi,FCTranSuc,Chi,SwapSuc,Chi,ReptSuc)
legend('Crank Shaft','Slide','Pivot','Full Chain Rot',...
        'Full Chain Trans','Swap','Reptate')
xlabel('Chi')
ylabel('Acceptance Rate')

figure(2)
hold on
plot(Chi,CrankShaftAmp,Chi,FCRotAmp,Chi,FCTranAmp)
legend('Crank Shaft','Full Chain Rog','Full Chain Trans');

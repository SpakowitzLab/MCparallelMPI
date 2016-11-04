function bindingCurve



dir='../data/';
%dir='../../C_0';
%dir='../../C_28/';
disp('reading params')
params=dlmread(strcat(dir,'paramsv1'));

NT=params(5)
LBox=params(10)
EU = params(11);
EM = params(12);
HP1_Bind = params(13);
LK = params(14);

RR=LBox/2;
out=dlmread(strcat(dir,sprintf('/out1v%d',1)),'',1,0);
snapMax=out(end,1);

snapMax=66;

reps=1:31;
nreps=length(reps);
snaps=50:snapMax;
nsnaps=length(snaps);


MU=zeros(nreps,1)*NaN;
MBcurve=zeros(nreps,1);
MUcurve=zeros(nreps,1);
WBcurve=zeros(nreps,1);
WUcurve=zeros(nreps,1);



for rep=reps
    out=dlmread(strcat(dir,sprintf('/out1v%d',rep)),'',1,0);
    MU(rep)=out(end,14);
    for snap=snaps
         
        r=dlmread(strcat(dir,sprintf('r%dv%d',snap,rep)));
        METH=logical(r(:,5));
        Bound=logical(r(:,4));

        %fprintf('rep=%d, snap=%d, nBound=%d, nMeth=%f\n',...
        %        rep,snap,sum(Bound),sum(METH));
      
        MBcurve(rep)=MBcurve(rep)+sum(METH & Bound);
        MUcurve(rep)=MUcurve(rep)+sum(METH & ~Bound);
        WBcurve(rep)=WBcurve(rep)+sum(~METH & Bound);
        WUcurve(rep)=WUcurve(rep)+sum(~METH & ~Bound);      

        if ((snap==snapMax(1,1)) & (rep==1))
                NMeth=sum(METH);
                NWeth=sum(~METH);
        end
    end
end
MBcurve=MBcurve/(nsnaps*NMeth); 
MUcurve=MUcurve/(nsnaps*NMeth);
WBcurve=WBcurve/(nsnaps*NWeth);
WUcurve=WUcurve/(nsnaps*NWeth);

figure(1)


plot(MU,MBcurve,'r'); hold on
plot(MU,WBcurve,'b');


mupts=MU(1):0.1:MU(end);
plot(mupts,1./(1+exp(-EM-mupts)),':r')
plot(mupts,1./(1+exp(-EU-mupts)),':b')
xlabel('\mu')
ylabel('P(bound)')
xlim([MU(1) MU(end)])
ylim([0 1])
legend('Meth','Un-Meth','Meth, no-int','Un-Meth, no-int')

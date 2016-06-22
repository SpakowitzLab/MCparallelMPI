%function [dum]=scalcfun(testnum,chemcp,lk)

testnum=1;
chemcp=2;
lk=50;

%%%%%%%%% PARAMS %%%%%%%%%%%%
boxl=15;dk=0.2;Ree=2.0;
FA=0.5;EPS=0.01;LAM=-0.75;G=5;
%lk=20;  %number of k calculations
lksample=20;
%%%%%%%%% END OF PARAMS %%%%%%

%%%%%%%%% INPUT %%%%%%%%%%
%chilistname=sprintf('/tower12/home/shifan/polymem-wlc-pt/rand-pt-03-06-15-%d-%d/',testnum,chemcp);
%SPINODAL=load('/tower12/home/shifan/polymem-wlc-pt/scalcbatch/analytical/chivals');
%CHIV=load([chilistname,'chilist'])
%CHIV=CHIV(1);

%NCHI=length(CHIV);
CHI=0;
NCHI=1;

snap0=1;snapf=51;skip=1;
snapv=snap0:skip:snapf;
nsnap=length(snapv);
%%%%%%%%% INPUT END %%%%%%%%%%

%%%% generate vectors on unit sphere %%%%
[kb,numk]=basisgen(lk,lksample,boxl);nk=length(kb);

%%%% start calculations %%%%
for cnum=1:NCHI
%CHI=CHIV(cnum);
col=(cnum-1)/(NCHI-1);

test=[];knum=1;
k=[];
%S=[];
S=zeros(nsnap,42);
Scntind=zeros(nsnap,1);
Stot=zeros(nk,1);

for kk=snapv
%for kk=1
disp(['chi=',num2str(CHI),' // snap#', num2str(kk)])
%cmat=rankpt(testnum,chemcp,kk);
%ind=find(abs(cmat(:,2)-CHI)<1e-2);
%index=cmat(ind,1);
%indexname=sprintf('/tower12/home/shifan/polymem-wlc-pt/rand-pt-03-06-15-%d-%d/rand-wlc-%d/',testnum,chemcp,index);
indexname='../';

    r=load([indexname,'data/r',num2str(kk)]);
    n=length(r);

    xr=r(:,1);yr=r(:,2);zr=r(:,3);
    xr=xr-boxl*floor(xr./boxl);
    yr=yr-boxl*floor(yr./boxl);
    zr=zr-boxl*floor(zr./boxl);

    id=r(:,4);
    f=sum(r(:,4)/n);

    t=2*(id-0.5);
    t=0.5*(t+1-2*f);

    jnum=1;
    for jj=snapv
    disp(['  chi=',num2str(CHI),' // snap#', num2str(jj)])
%    cmat=rankpt(testnum,chemcp,jj);
%    ind=find(abs(cmat(:,2)-CHI)<1e-2);
%    index=cmat(ind,1);
%    indexname=sprintf('/tower12/home/shifan/polymem-wlc-pt/rand-pt-03-06-15-%d-%d/rand-wlc-%d/',testnum,chemcp,index);
indexname='../';

        r=load([indexname,'data/r',num2str(jj)]);
        n=length(r);

        xr=r(:,1);yr=r(:,2);zr=r(:,3);
        xr2=xr-boxl*floor(xr./boxl);
        yr2=yr-boxl*floor(yr./boxl);
        zr2=zr-boxl*floor(zr./boxl);

        id=r(:,4);
        f=sum(r(:,4)/n);

        t=2*(id-0.5);
        t2=0.5*(t+1-2*f);

    for ii=1:nk
%      ii
      xcom=kb(ii,1)*xr;
      ycom=kb(ii,2)*yr;
      zcom=kb(ii,3)*zr;

      xcom2=kb(ii,1)*xr2;
      ycom2=kb(ii,2)*yr2;
      zcom2=kb(ii,3)*zr2;

      Scos=(cos(xcom+ycom+zcom).*t);
      Ssin=(sin(xcom+ycom+zcom).*t);

      Scos2=(cos(xcom2+ycom2+zcom2).*t2);
      Ssin2=(sin(xcom2+ycom2+zcom2).*t2);

%      Stot(ii)=sum(Scos).^2+sum(Ssin).^2;
      Stot(ii)=sum(Scos).*sum(Scos2)+sum(Ssin).*sum(Ssin2);
    end

%    [k,S(knum,:)]=scalcavg(kb,Stot);
    [k,Stemp]=scalcavg(kb,Stot);
    indS=abs(knum-jnum)+1;
    S(indS,:)=S(indS,:)+transpose(Stemp);
    Scntind(indS)=Scntind(indS)+1;

    jnum=jnum+1; %each snapshot
    end

    knum=knum+1;
end

S=S./n;  %normalization

for nsnap=1:nsnap
  S(nsnap,:)=S(nsnap,:)/Scntind(nsnap);
end

%Sa=mean(S,1);

%whos k Sa
%test=[k*Ree,transpose(Sa)];
%filename=sprintf('SMC_SIM%dCHEM%dCHI%.8f',testnum,chemcp,CHI);
%filename=sprintf('SMC_SIM%dCHEM%dCHI%.8f',testnum,chemcp,CHI);
%dlmwrite(filename,test,'delimiter','\t','precision',3)

sizeS=size(S);
nk=sizeS(2);

%for nkmode=1:nk
%  S(:,nkmode)=S(:,nkmode)./S(1,nkmode);
%end

figure;hold
for nkmode=1:5:nk
  col = (nkmode-1)/(nk-1);
plot(S(:,nkmode),'color',[col 0 1-col],'linewidth',2);
end

figure;hold
for ii=1:nsnap
  ii
  col = (ii-1)/(nsnap-1);
  plot(k*Ree,S(ii,:),'color',[col 0 1-col])
end

end

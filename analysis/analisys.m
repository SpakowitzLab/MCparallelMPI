function analisys

    %dir='../firstInrData5_30_16/';
    dir='../data/';
    disp('reading params')
    params=dlmread(strcat(dir,'paramsv1')); % parameters
    %NT=params(1);
    %N=params(2);
    %NB=params(3);
    %NP=params(4);
    NT=params(5);
    %G=params(6);
    %L0=params(7);
    %CHI0=params(8);
    %Fpoly=params(9);
    LBox=params(10);
    EU = params(11);
    EM = params(12);
    HP1_Bind = params(13);
    LK = params(14);
    %FRMMETH = params(15);
    %F_METH = params(16);
    %LAM_METH = params(17);
    R=LBox/2;
    
mu_data=dlmread(strcat(dir,'mu'));
nrep=length(mu_data);

nbins=50;
edges=[0:(nbins)].*1/(nbins);
densitymb=zeros(nbins,nrep);
densitymu=zeros(nbins,nrep);
densitywb=zeros(nbins,nrep);
densitywu=zeros(nbins,nrep);



fmb=zeros(1,nrep);
fmu=zeros(1,nrep);
fwb=zeros(1,nrep);
fwu=zeros(1,nrep);


Iseq=115:1:120;
naverage=length(Iseq);
for rep=1:nrep
    for I=Iseq
        string=strcat(dir,sprintf('u0v%d',rep));
        u0=dlmread(string); % u0
        METH=logical(u0(:,5));
        string=strcat(dir,sprintf('r%dv%d',I,rep));
        
        data=dlmread(string);
        bound=logical(data(:,4));
        
        mb=data(  (METH.*   bound)==1,1:3);
        mu=data(  (METH.* (~bound))==1,1:3);
        wb=data(((~METH).*  bound)==1,1:3);
        wu=data(((~METH).*(~bound))==1,1:3);        
        
            
        figure(1)
        fracV=(((mb(:,1)-R).^2 + (mb(:,2)-R).^2  + (mb(:,3)-R).^2).^(3/2))./R^3;
        h=histogram(fracV,edges);
        densitymb(:,rep)=densitymb(:,rep)+h.Values';
        fmb(1,rep)=fmb(1,rep)+sum(h.Values);
        fmbFluk(I,rep)=sum(h.Values);

        fracV=(((mu(:,1)-R).^2 + (mu(:,2)-R).^2  + (mu(:,3)-R).^2).^(3/2))./R^3;
        h=histogram(fracV,edges);
        densitymu(:,rep)=densitymu(:,rep)+h.Values';
        fmu(1,rep)=fmu(1,rep)+sum(h.Values);
        fmuFluk(I,rep)=sum(h.Values);
        

        fracV=(((wb(:,1)-R).^2 + (wb(:,2)-R).^2  + (wb(:,3)-R).^2).^(3/2))./R^3;
        h=histogram(fracV,edges);
        densitywb(:,rep)=densitywb(:,rep)+h.Values';
        fwb(1,rep)=fwb(1,rep)+sum(h.Values);
        fwbFluk(I,rep)=sum(h.Values);
        

        fracV=(((wu(:,1)-R).^2 + (wu(:,2)-R).^2  + (wu(:,3)-R).^2).^(3/2))./R^3;
        h=histogram(fracV,edges);
        densitywu(:,rep)=densitywu(:,rep)+h.Values';
        fwu(1,rep)=fwu(1,rep)+sum(h.Values);
        fwuFluk(I,rep)=sum(h.Values);
        

        
    end
        
end

fmb=fmb*nrep/(NT*naverage);
fmu=fmu*nrep/(NT*naverage);
fwb=fwb*nrep/(NT*naverage);
fwu=fwu*nrep/(NT*naverage);

figure(2)
mupts=-3:0.1:4;
ftheory=exp( -(-EM-mupts))./(1+exp( -(-EM-mupts)));
plot(mupts,ftheory,'-k',mu_data,fmb./(fmb+fmu),'-b')
xlabel('mu')
ylabel('fraction bound')
hold on
ftheory=exp( -(-EU-mupts))./(1+exp( -(-EU-mupts)));
plot(mupts,ftheory,'-k',mu_data,fwb./(fwb+fwu),'-g')
legend('meth no intr.','meth avj','un-meth no intr','un-meth avj')
title('binding curve')

for I=Iseq
    plot(mu_data,fmbFluk(I,:)./(fmbFluk(I,:)+fmuFluk(I,:)),'.b')
    plot(mu_data,fwbFluk(I,:)./(fwbFluk(I,:)+fwuFluk(I,:)),'.g')
end


densityTotal=densitywu+densitywb+densitymu+densitymb;

repPick=[1,8,11,12,15,39];

for p=1:6
    xpts=0.5*(edges(1:end-1)+edges(2:end));
    figure(3)
    factor=nbins/(NT*naverage*4*pi()/3);
    subplot(2,3,p)
    hold on
    rep=repPick(p);
    plot(xpts.^0.33333,densityTotal(:,rep)*factor,'-b')
    plot(xpts.^0.33333,densitywu(:,rep)*factor,'-om')
    plot(xpts.^0.33333,densitywb(:,rep)*factor,'-og')
    plot(xpts.^0.33333,densitymu(:,rep)*factor,'-sr')
    plot(xpts.^0.33333,densitymb(:,rep)*factor,'-sc')
    Rvec=(0:200)/200;
    Yvec=(sin(pi()*Rvec).^2)./(2*pi()*Rvec.^2);
    plot(Rvec,Yvec,'-k')
    string=sprintf('mu=%f',mu_data(rep));
    title(string)
    xlabel('r/R')
    ylabel('density')
    set(gca,'YTickLabel','')
    if p==1
        legend('totoal','un-meth, un-bound','un-meth, bound', ...
              'meth, un-bound','meth-bound', 'gaussian')
    end
    
end
    %string=sprintf('Total density. LK=%f, R=%f',LK,R);
    %title(string)

end
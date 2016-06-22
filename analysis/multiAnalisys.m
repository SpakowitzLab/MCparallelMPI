function multiAnalisys

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

for set=1:7
    Iseq=105+10*set:1:110+10*set;
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

            fracV=(((mu(:,1)-R).^2 + (mu(:,2)-R).^2  + (mu(:,3)-R).^2).^(3/2))./R^3;
            h=histogram(fracV,edges);
            densitymu(:,rep)=densitymu(:,rep)+h.Values';
            fmu(1,rep)=fmu(1,rep)+sum(h.Values);


            fracV=(((wb(:,1)-R).^2 + (wb(:,2)-R).^2  + (wb(:,3)-R).^2).^(3/2))./R^3;
            h=histogram(fracV,edges);
            densitywb(:,rep)=densitywb(:,rep)+h.Values';
            fwb(1,rep)=fwb(1,rep)+sum(h.Values);


            fracV=(((wu(:,1)-R).^2 + (wu(:,2)-R).^2  + (wu(:,3)-R).^2).^(3/2))./R^3;
            h=histogram(fracV,edges);
            densitywu(:,rep)=densitywu(:,rep)+h.Values';
            fwu(1,rep)=fwu(1,rep)+sum(h.Values);



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

end

end
function SpaceCorrilation
    dir='../firstInrData5_30_16/';
   disp('reading params')
    params=dlmread(strcat(dir,'paramsv1')); % parameters
    LBox=params(10);
    NT=params(1);
    R=LBox/2;


mu_data=dlmread(strcat(dir,'mu'));
nrep=length(mu_data);


% for corrilation plot
nbins2=100;
edges2=(0:(nbins2)).*1/(nbins2);
edges2=edges2;
MMr3hist=zeros(nbins2,nrep);
WWr3hist=zeros(nbins2,nrep);
MWr3hist=zeros(nbins2,nrep);

skip=20;
repSkip=5;

naverage=0;
for rep=1:repSkip:nrep
    for I=50:5:95
        fprintf('No int. rep=%d, I=%d\n',rep,I)
        naverage=naverage+1;
        string=sprintf(strcat(dir,'u0v%d'),rep);
        u0=dlmread(string); % u0
        METH=logical(u0(:,5));
        string=sprintf(strcat(dir,'r%dv%d'),I,rep);
        data=dlmread(string);
        bound=logical(data(:,4));
        
        mb=data(  (METH.*   bound)==1,1:3);
        mu=data(  (METH.* (~bound))==1,1:3);
        wb=data(((~METH).*  bound)==1,1:3);
        wu=data(((~METH).*(~bound))==1,1:3);        
        
        if (size(mb,1)+size(mu,1)+size(wb,1)+size(wu,1))~=NT
            error('NT is wrong')
        end  
        figure(1)
        

        mth=[mb;mu];
        wth=[wb;wu];
        Nmth=size(mth,1);
        Nwth=size(wth,1);
        % Calculate MM corrilation
        MMcor=zeros(ceil((Nmth/skip)*((Nmth/skip)-1)/2),1);
        kk=1;
        for ii=1:skip:Nmth
           for jj=(ii+1):skip:Nmth
               MMcor(kk,rep)=norm(mth(ii,:)-mth(jj,:));
               kk=kk+1;
           end
        end
        h=histogram((MMcor(1:kk-1,rep)./R).^3,edges2);
        MMr3hist(:,rep)=MMr3hist(:,rep)+h.Values';

        % Calculate WW corrilation
        WWcor=zeros(ceil((Nwth/skip)*((Nwth/skip)-1)/2),1);
        kk=1;
        for ii=1:skip:Nwth
           for jj=(ii+1):skip:Nwth
               WWcor(kk,rep)=norm(wth(ii,:)-wth(jj,:));
               kk=kk+1;
           end
        end
        h=histogram((WWcor(1:kk-1,rep)./R).^3,edges2);
        WWr3hist(:,rep)=WWr3hist(:,rep)+h.Values';

        % Calculate MWcorrilation
        MWcor=zeros(ceil(Nwth*Nmth/(skip^2)));
        kk=1;
        for ii=1:skip:Nwth
           for jj=1:skip:Nmth
               MWcor(kk,rep)=norm(wth(ii,:)-mth(jj,:));
               kk=kk+1;
           end
        end
        h=histogram((MWcor(1:kk-1,rep)./R).^3,edges2);
        MWr3hist(:,rep)=MWr3hist(:,rep)+h.Values';

        
    end
        
end

fprintf('number M/number W=%f\n',Nmth/Nwth)
fprintf('number M^2/number W^2=%f\n',Nmth^2/Nwth^2)

for rep=1:repSkip:nrep
    figure(3)
    xpts=0.5*(edges2(1:end-1)+edges2(2:end));
    h1=plot(xpts*R,MMr3hist(:,rep),'k');
    hold on
    title('Spatial corrilation')
end

for rep=1:repSkip:nrep
    figure(3)
    xpts=0.5*(edges2(1:end-1)+edges2(2:end));
    h2=plot(xpts*R,WWr3hist(:,rep),'r');
    hold on
    title('Spatial corrilation')
end

for rep=1:repSkip:nrep
    figure(3)
    xpts=0.5*(edges2(1:end-1)+edges2(2:end));
    h3=plot(xpts*R,MWr3hist(:,rep),'y');
    hold on
    title('Spatial corrilation')
end
% ---------------------------------------------------------
%
%    After turning on interaction
%
% -------------------------------------------------------
MMr3hist=zeros(nbins2,nrep);
WWr3hist=zeros(nbins2,nrep);
MWr3hist=zeros(nbins2,nrep);
naverage=0;
for rep=1:repSkip:nrep
    for I=150:5:195
        fprintf('Int. on rep=%d, I=%d\n',rep,I)
        naverage=naverage+1;
        string=strcat(dir,sprintf('u0v%d',rep));
        u0=dlmread(string); % u0
        METH=logical(u0(:,5));
        string=sprintf(strcat(dir,'r%dv%d'),I,rep);
        data=dlmread(string);
        bound=logical(data(:,4));
        
        mb=data(  (METH.*   bound)==1,1:3);
        mu=data(  (METH.* (~bound))==1,1:3);
        wb=data(((~METH).*  bound)==1,1:3);
        wu=data(((~METH).*(~bound))==1,1:3);        
        
        if (size(mb,1)+size(mu,1)+size(wb,1)+size(wu,1))~=NT
            error('NT is wrong')
        end  
        figure(1)
        
        mth=[mb;mu];
        wth=[wb;wu];
        Nmth=size(mth,1);
        Nwth=size(wth,1);
        % Calculate MM corrilation
        MMcor=zeros(ceil((Nmth/skip)*((Nmth/skip)-1)/2),1);
        kk=1;
        for ii=1:skip:Nmth
           for jj=(ii+1):skip:Nmth
               MMcor(kk,rep)=norm(mth(ii,:)-mth(jj,:));
               kk=kk+1;
           end
        end
        h=histogram((MMcor(1:kk-1,rep)./R).^3,edges2);
        MMr3hist(:,rep)=MMr3hist(:,rep)+h.Values';

        % Calculate WW corrilation
        WWcor=zeros(ceil((Nwth/skip)*((Nwth/skip)-1)/2),1);
        kk=1;
        for ii=1:skip:Nwth
           for jj=(ii+1):skip:Nwth
               WWcor(kk,rep)=norm(wth(ii,:)-wth(jj,:));
               kk=kk+1;
           end
        end
        h=histogram((WWcor(1:kk-1,rep)./R).^3,edges2);
        WWr3hist(:,rep)=WWr3hist(:,rep)+h.Values';

        % Calculate MWcorrilation
        MWcor=zeros(ceil(Nwth*Nmth/(skip^2)));
        kk=1;
        for ii=1:skip:Nwth
           for jj=1:skip:Nmth
               MWcor(kk,rep)=norm(wth(ii,:)-mth(jj,:));
               kk=kk+1;
           end
        end
        h=histogram((MWcor(1:kk-1,rep)./R).^3,edges2);
        MWr3hist(:,rep)=MWr3hist(:,rep)+h.Values';

        
    end
        
end

fprintf('number M/number W=%f\n',Nmth/Nwth)
fprintf('number M^2/number W^2=%f\n',Nmth^2/Nwth^2)

for rep=1:repSkip:nrep
    figure(3)
    xpts=0.5*(edges2(1:end-1)+edges2(2:end));
    h4=plot(xpts*R,MMr3hist(:,rep),'g');
    hold on
    title('Spatial corrilation')
end

for rep=1:repSkip:nrep
    figure(3)
    xpts=0.5*(edges2(1:end-1)+edges2(2:end));
    h5=plot(xpts*R,WWr3hist(:,rep),'b');
    hold on
    title('Spatial corrilation')
end
for rep=1:repSkip:nrep
    figure(3)
    xpts=0.5*(edges2(1:end-1)+edges2(2:end));
    h6=plot(xpts*R,MWr3hist(:,rep),'c');
    hold on
    title('Number of neighbors')
end
list={'Int off MM','Int off WW','Int off MW', ...
      'Int on MM','Int on WW','Int on MW'};
legend([h1,h2,h3,h4,h5,h6],list)
set(gca,'YTickLabel',[])

end
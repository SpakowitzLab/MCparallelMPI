function [k,S, varargout]=scalc(r,boxl,lksample)
nout=max(nargout,1);
%% calculate structure factor]
% INPUTS: 
% OUTPUT: 

%% %%%%%%%%% PRE-PROCESS DATA %%%%%%%%%
n=length(r);      % number of beads
id=r(:,4);        % chemical identities
f=sum(r(:,4)/n);  % fraction of A-type monomers

%% apply periodic boundary conditions to coordinates
xr=r(:,1);yr=r(:,2);zr=r(:,3);
xr=xr-boxl*floor(xr./boxl);
yr=yr-boxl*floor(yr./boxl);
zr=zr-boxl*floor(zr./boxl);

%% transform chemical identity to \sigma_A-f_A
t=2*(id-0.5);
t=0.5*(t+1-2*f);

%%%% start calculations %%%%
ncomplete=8;
kcomplete=(2*pi/boxl)*ncomplete;
[kb]=basisgen(lksample,boxl,ncomplete);  % create wavevectors
nk=length(kb);               % number of wavevectors
Stot=zeros(nk,1);
for ii=1:nk
  xcom=kb(ii,1)*xr;
  ycom=kb(ii,2)*yr;
  zcom=kb(ii,3)*zr;

  Scos=(cos(xcom+ycom+zcom).*t);
  Ssin=(sin(xcom+ycom+zcom).*t);
  Stot(ii)=sum(Scos).^2+sum(Ssin).^2;
end

%%%%%%%%% POST-PROCESS DATA %%%%%%%%%
[k,S]=scalcavg(kb,Stot);   % rotational average S
S=S./n;                    % include bead normalization

if nout==2
    return
elseif nout==3
    varargout{1}=kcomplete;
    return    
elseif nout==5
    varargout{1}=kcomplete;
    varargout{2}=kb;
    varargout{3}=Stot./n;
    return
else
    error('wrong number of return values')
end

end
%% generate basis wavevectors
function [kb]=basisgen(Nsample,L,ncomplete)
%Ksample = 200;
rng(123,'twister');
test = unifrnd(0,1,(Nsample+1)^3,1);
rng('shuffle')
%test=ones((Nsample+1)^3,1);
kb=zeros((Nsample+1)^3,4);
ii=0;
jj=0;
for n1=0:Nsample
  for n2=0:Nsample
    for n3=0:Nsample
        ii=ii+1;
        if (n1+n2+n3)==0
            continue
        end
        p=(n1*n1+n2*n2+n3*n3)/(ncomplete*ncomplete);
        p=(1/p)*exp(-sqrt(p) +1);
        if p < test(ii)
            continue
        end
        jj=jj+1;
        k1=(2*pi/L)*n1;
        k2=(2*pi/L)*n2;
        k3=(2*pi/L)*n3;
        km=sqrt(k1*k1+k2*k2+k3*k3);
        kb(jj,:)=[k1,k2,k3,km];
    end
  end
end
kb=kb(1:jj,:);


%only select some k values at high k
% nk=length(unique(kb(:,4)));
% kmag=unique(kb(:,4));
% ind =unique(round(logspace(0,log10(nk),Ksample)));
% 
% kbnew=[];
% for ii=1:length(ind)
%   kbnew=[kbnew;kb(kb(:,4)==kmag(ind(ii)),:)];
% end

end
%% average structure factor
function [k,S]=scalcavg(kb,Stot)


% tol=1e-3;
% k=round(kb/tol)*tol;
% k=unique(k(:,4));
% 
% numk=length(k);
% S=zeros(numk,1);
% 
% for ii=1:numk
%    index=find(abs(kb(:,4)-k(ii))<=tol);
%    S(ii)=mean(Stot(index));
%  %    S(ii)=mean(Stot(abs((kb(:,4)-k(ii))./k(ii))<=tol));
% end
% sprintf('max(km)=%f in average',max(k))


%sort kb and Stot
tol=10^-3;
[ks,I]=sort(kb(:,4));
kb=kb(I,:);
Stot=Stot(I);

n=0;
avj=0;
nk=length(kb);
k=zeros(nk,4)*NaN;
S=zeros(nk,1)*NaN;
jj=0;
for ii=1:nk
    if n==0
        n=1;
        avj=Stot(ii);
    elseif ((ks(ii)-ks(ii-1))/ks(ii) <tol)
        n=n+1;
        avj=avj+Stot(ii);
    else
        jj=jj+1;
        k(jj,:)=ks(ii-1);
        S(jj)=avj/n;
        
        n=1;
        avj=Stot(ii);
    end
end
jj=jj+1;
k(jj,:)=ks(ii-1);
S(jj)=avj/n;
S=S(1:jj);
k=k(1:jj,4);

end

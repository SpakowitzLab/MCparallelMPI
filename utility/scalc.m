function [k,S]=scalc(r,boxl,lksample)
%% calculate structure factor]
% INPUTS: 
% OUTPUT: 
kb=basisgen(lksample,boxl);

%%%%%%%%% PRE-PROCESS DATA %%%%%%%%%
n=length(r);      % number of beads
nk=length(kb);    % number of wavevectors
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

end

function [kb]=basisgen(Nsample,L)
% generate basis wavevectors

kb=[];
for n1=0:Nsample
  for n2=0:Nsample
    for n3=0:Nsample
      if (n1+n2+n3)>0
        k1=(2*pi/L)*n1;
        k2=(2*pi/L)*n2;
        k3=(2*pi/L)*n3;
        km=sqrt(k1*k1+k2*k2+k3*k3);
        kb=[kb;k1,k2,k3,km];
      end
    end
  end
end
end

function [k,S]=scalcavg(kb,Stot)
% average structure factor

tol=1e-3;
k=round(kb/tol)*tol;
k=unique(k(:,4));

numk=length(k);
S=zeros(numk,1);

for ii=1:numk
%   index=find(abs(kb(:,4)-k(ii))<=tol);
%   S(ii)=mean(Stot(index));
    S(ii)=mean(Stot(abs(kb(:,4)-k(ii))<=tol));
end
end
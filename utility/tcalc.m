function [cosTheta,pdf,n,S]=tcalc(r,NP,L0)

%% calculate end-to-end distribution of polymers

%% find end-end vectors
NB=length(r);
N=NB/NP;
L=(N-1)*L0;
I1=((1:1:NP)-1)*N+1;
I2=((1:1:NP)-1)*N+N;
REND=(r(I2,1:3)-r(I1,1:3))./L;

%% find sphericity
num=0;denom=0;
for ii=1:size(REND,1)
    num=num+REND(ii,:)'*REND(ii,:);
    denom=denom+REND(ii,:)*REND(ii,:)';
end
M=num./denom;
[V,D]=eig(M);

vec=[D(1,1); D(2,2); D(3,3)];
[vec,ind]=sort(vec,'descend');


S=(3/2)*(vec(2)+vec(3)); % Sphericity
n=V(:,ind(1)); % Primary axis

mag=sqrt(sum(REND.^2,2)); % Magnitude of end-end distances

co=abs(REND*n)./mag; % cos with primary axis

edges=0:0.05:1;
[pdf,bins]=histcounts(co,edges,'Normalization','pdf'); 
cosTheta=0.5*(bins(2:end)+bins(1:end-1));

end
path=sprintf('../../data/');
npts=100000;
nrep=7;
chiValues=zeros(1,nrep);
xchi=zeros(npts,nrep); % # save points x # terms
%chiValues=[0.0143,0.3215,10.313,20.979,30.313,40.312,50.313,60];
for rep=1:nrep
     data=dlmread(strcat(path,sprintf('rv%d',rep)));
     xchi(:,rep)=data(1:npts,1);
     data=dlmread(strcat(path,sprintf('cofv%d',rep)));
     chiValues(1,rep)=data(1,1);
end
for rep=1:nrep
    figure(rep)
    nbin=30;
    [N,edges]=histcounts(xchi(:,rep),nbin);
    centers=0.5*(edges(1:end-1)+edges(2:end));
    N=N/(npts*(centers(2)-centers(1)));
    hold on
    theory=exp(-chiValues(rep)*centers.^2)*sqrt(chiValues(rep)/pi);
    plot(centers,N,'o',centers,theory)
end
disp(chiValues)
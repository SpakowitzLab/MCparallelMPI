function radialDensity

dir='../data/';
%dir='../../C_0/';
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

fprintf('NT=%d, LBox=%f, C=%f, L_khun=%f\n',NT,LBox,HP1_Bind,LK)
% out=dlmread(strcat(dir{condition},sprintf('/out1v%d',rep)),'',1,0);

nbins=40;
edges=(0:nbins)*RR/nbins;

Vol=(4/3)*pi()*(edges(2:end).^3)-(4/3)*pi()*(edges(1:end-1).^3);

centers=0.5*(edges(2:end)+edges(1:end-1));

reps=1:31;
nreps=length(reps);
snaps=40:52;
nsnaps=length(snaps);

radialDist=zeros(nreps,nbins);

for rep=reps
    fprintf('rep %d of %d\n',rep,nreps)
    for snap=snaps
        r=dlmread(strcat(dir,sprintf('r%dv%d',snap,rep)));

        R=r(:,1:3);
        R=R-ones(NT,3)*(LBox/2);
        R=sqrt(sum(R.^2,2));
        %fprintf('snap=%d, rep=%d, mean(R)=%f, max(R)=%f\n',snap,rep,mean(R),max(R));
        counts=histcounts(R,edges);
        radialDist(rep,:)=radialDist(rep,:)+counts;
    end
end


figure (1) % Density distribution
for rep=reps
        col=rep/nreps;

        %density=radialDist(snap,:)./Vol;
        %density=density/nreps;
        prob=(nbins/RR)*radialDist(rep,:)/sum(radialDist(rep,:));
        figure(1)
        plot(centers,prob,'color',[col 0 1-col]); hold on
        figure(2)
        plot(centers,prob*NT*0.1./(4*pi()*centers.^2),'color',[col 0 1-col]); hold on
end
figure(1)
rpts=(1:100)/100;
plot(RR*rpts,2*(sin(pi()*rpts)).^2/RR,'k-')

figure(2)
rpts=(1:100)/100;
plot(RR*rpts,(2.*(sin(pi()*rpts)).^2/RR).*(NT*0.1./(4*pi()*RR*RR*rpts.^2)),'k-')
ylim([0 1])
xlabel('r/R')
ylabel('density')


end

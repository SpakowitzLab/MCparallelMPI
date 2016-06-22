function pplot(N,EPS)
dir = 'pdata/';
FNUM=round(100*EPS);

XG=transpose(linspace(0,1,5001));
if FNUM<=2000  % use analytical wormlike chain statistics
    G=load([dir,'out',int2str(FNUM),'.txt']);
    plot(XG,G.*power(XG,2)*(4*pi),'k-','LineWidth',2)
else  % use Gaussian chain statistics
    DXG=XG(2)-XG(1);
    G=exp(-1.5*N*XG.^2).*XG.*XG;
    G=G./(sum(G).*DXG);
    plot(XG,G,'k-','LineWidth',2)
end

axis([0 1 0 10])
set(gca,'FontSize',14)
xlabel('R/L','FontSize',18)
ylabel('P(R)','FontSize',18)
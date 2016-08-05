function DirectionPlot

path=sprintf('../../simData/lamNeg075_r7_27_16/');
boxl = 20;
lksample = 20;
NP = 2000;
EPS = 0.01;
G = 5;
N = 8;
L0=0.40667; %2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5)
Ree = 2;


repMIN=20;
repSkip=1;
repMAX=55;
snapMin=26;
snapSkip=4;
snapMax=75;

nrep=length(repMIN:repSkip:repMAX);
v=0;
chiVec={};

TEST1=1;
TEST5=0;
TEST6=1;
TEST7=1;

for rep=repMIN:repSkip:repMAX
    v=v+1;
    first=1;
    navj=0;
    fprintf('calculating replica %d of %d ...\n',v,nrep)
    for snap=snapMin:snapSkip:snapMax
        r=dlmread(strcat(path,sprintf('r%dv%d',snap,rep)));
        if (TEST1)
            % end-to-end distribution
            if first
                [X,P]=pcalc(r,NP,L0);
            else
                [dX,dP]=pcalc(r,NP,L0);
                X=X+dX;
                P=P+dP;
            end
        end
        if (TEST6)
            % end-to-end distribution
            if first
                [cosTheta,pdf,nvec,S]=tcalc(r,NP,L0);
            else
                [dcosTheta,dpdf,dnvec,dS]=tcalc(r,NP,L0);
                if any(dcosTheta~=cosTheta)
                    error('Error in Direction Plot 6')
                end
                pdf=pdf+dpdf;
                S=S+dS;
            end
        end
        if (TEST5)
            if first
                [k,s,kcomplete,kdir,sdir]=scalc(r,boxl,lksample);
            else
                [dk,ds,kcomplete,dkdir,dsdir]=scalc(r,boxl,lksample);
                if any(dk~=k)
                    error('Error in Direction Plot 5')
                end
                s=s+ds;
                kdir=kdir+dkdir;
                sdir=sdir+dsdir;
            end
        end
        first=0;
        navj=navj+1;
    end
    if (repMAX-repMIN>1)
        col=(rep-repMIN)/(repMAX-repMIN);
    else
        col=1;
    end
    out1=dlmread(strcat(path,sprintf('out1v%d',rep)),'',1,0);
    chiVec{v}=num2str([out1(end,13), rep]);
%% PLOT
    if (TEST1)
        X=X/navj;
        P=P/navj;
        figure(1);
        plot(X,P,'color',[col 0 1-col]);hold on;
    end
    if(TEST5)
        s=s/navj;
        figure(5);
        plot(k*Ree,s,'color',[col 0 1-col]);hold on;
    end
    if(TEST6)
        pdf=pdf/navj;
        figure(6);
        plot(cosTheta,pdf,'color',[col 0 1-col]);hold on;
    end
    if(TEST7)
        figure(7)
        S=S/navj;
        plot(out1(end,13),S,'o','color',[col 0 1-col]);hold on;
    end
%     size(sdir)
%     size(kdir)
%     [peakH,peakI]=max(sdir);
    
end

%% PLOT Theory, legends, and axis
if (TEST1)
    figure(1); hold on;
    pplot(N*G,EPS*N*G)
    xlabel('R/L');ylabel('P(R/L)')
    legend(chiVec)
end
if (TEST5)
    figure(5); hold on;
%     filename = sprintf('sdata/Seps%.3flam%.2f',EPS,LAM);
%     S = load(filename);
%     KS = S(:,1);
%     SS = 1./(-2*CHI+1./S(:,2));
%     plot(KS,SS);
    xlabel('q*Ree');ylabel('S(q)')
    plot(kcomplete*Ree,0.01,'o');
    set(gca,'xscale','log');set(gca,'yscale','log')
    legend(chiVec)
end
if (TEST6)
    figure(6); hold on;
    xlabel('cos(\theta)'); ylabel('pdf');
    legend(chiVec)
end
if (TEST7)
    figure(7);hold on;
    xlabel('\chi'); ylabel('Sphericity');
end

path=sprintf('../MCPoly/MCparallelMPI/data');

%% parameters
boxl = 20;
DEL = 1;
V = 0.1;
lksample = 20;
NP = 2000;
EPS = 1.01;
G = 5;
N = 8;
L0=40667; %2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5)
LP=L0/(2*EPS);
Ree = 2;
LAM = -0.75;
FA = 0.5;
CHI = 0;
KAP = 10;

repMIN=1;
repSkip=5;
repMAX=7;
snapMin=1;
snapSkip=1;
snapMax=1;

%% TEST OPTIONS
TEST1 = 1;  % end-to-end distribution
TEST2 = 0;  % bead-bead distribution
TEST3 = 0;  % total and partial density
TEST4 = 0;  % radial distrbution function
TEST5 = 0;  % structure factor
TEST6 = 0;  % chemical potential (widom insertion)

%% Main Loop
%%
nrep=length(repMIN:repSkip:repMAX);
v=0;
chiVec={};
for rep=repMIN:repSkip:repMAX
    v=v+1;
    first=1;
    navj=0;
    fprintf('calculating %d of %d ...\n',v,nrep)
    for snap=snapMin:snapSkip:snapMax
        r=dlmread(strcat(path,sprintf('/r%dv%d',snap,rep)));
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
        if (TEST2)
            if first
                R = rcalc(r,NP);
            else
                R=R+ rcalc(r,NP);
            end
        end
        if (TEST3)
            if first
                [PHIA,PHIB] = r_to_phi(r,boxl,DEL,V);
            else
                [dPHIA,dPHIB] = r_to_phi(r,boxl,DEL,V);
                PHIA=PHIA+dPHIA;
                PHIB=PHIB+dPHIB;
            end
        end
        if (TEST5)
            if first
                [k,s]=scalc(r,boxl,lksample);
            else
                [dk,ds]=scalc(r,boxl,lksample);
                k=k+dk;
                s=s+ds;
            end
        end
        first=0;
        navj=navj+1;
    end
    col=(rep-repMIN)/(repMAX-repMIN);
    out1=dlmread(strcat(path,sprintf('/out1v%d',rep)));
    chiVec{v}=num2str(out1(end,11));
%% PLOT
    if (TEST1)
        X=X/navj;
        P=P/navj;
        figure(1);
        hold on;
        plot(X,P,'color',[col 0 1-col]);
    end
    if (TEST2)
        R=R/navj;
        NK = (0:N*G-1)'/G;
        figure(2); hold on
        plot(NK,sqrt(R)/Ree,'o','color',[col 0 1-col]);
    end
    if(TEST3)
        PHIA=PHIA/navj;
        PHIB=PHIB/navj;
        figure(3); hold on
        subplot(ceil(sqrt(nrep)),ceil(sqrt(nrep)),v);
        hist(PHIA+PHIB,50,'color',[col 0 1-col]);
        xlabel('\phi_A+\phi_B');ylabel('P(\phi)')
    end
    if(TEST5)
        k=k/navj;
        s=s/navj;
        figure(5);hold on
        plot(k*Ree,s,'color',[col 0 1-col]);
    end
end

%% PLOT Theory
if (TEST1)
    figure(1);
    hold on;
    pplot(N*G,EPS*N*G)
    xlabel('R/L');ylabel('P(R/L)')
    legend(chiVec)
end
if (TEST2)
    % bead-bead distribution
    R = rcalc(r,NP);
    NK = (0:N*G-1)'/G;
    NM = EPS*G;
    DM = EPS*NK*G;
    RM = NM-(1/2)*(1-exp(-2*NM));
    RD = DM-(1/2)*(1-exp(-2*DM));

    figure(2); hold on
    plot(NK,sqrt(RD/RM),'-')
    set(gca,'xscale','log');set(gca,'yscale','log')
    xlabel('J/G');ylabel('R_J/R_M')
end
if (TEST5)
    % load analytical theory
    figure(5); hold on
    filename = sprintf('sdata/Seps%.3flam%.2f',EPS,LAM);
    S = load(filename);
    KS = S(:,1);
    SS = 1./(-2*CHI+1./S(:,2));
    plot(KS,SS);
    xlabel('q');ylabel('S(q)')
    set(gca,'xscale','log');set(gca,'yscale','log')
    legend(chiVec)
end


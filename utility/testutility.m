clear;
close all
% inputs:
% r101, an example file for bead coordinates
% pdata, analytical solution for end-to-end distribution
%   pdata available shifan@tower12:~/testdir/pdata
% sdata, analytical solution of structure factors
%   sdata available shifan@tower12:~/testdir/sdata

%% read data
r = dlmread('../../data/r15v1');

%% TEST OPTIONS
TEST1 = 1;  % end-to-end distribution
TEST2 = 1;  % bead-bead distribution
TEST3 = 1;  % total and partial density
TEST4 = 1;  % radial distrbution function
TEST5 = 1;  % structure factor
TEST6 = 1;  % chemical potential (widom insertion)

%% parameters
boxl = 20;
DEL = 1;
V = 0.1;
lksample = 20;
NP = 2000;
EPS = 0.01;
G = 5;
N = 8;
L0=2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5);
LP=L0/(2*EPS);
Ree = 2;
LAM = -0.75;
FA = 0.5;
CHI = 0/G;
KAP = 10;

%% Single chain conformation tests
if (TEST1)
    % end-to-end distribution
    [X,P]=pcalc(r,NP,L0);
    figure;hold;
    plot(X,P);pplot(N*G,EPS*N*G)
    xlabel('R/L');ylabel('P(R/L)')
end

if (TEST2)
    % bead-bead distribution
    R = rcalc(r,NP);
    NK = (0:N*G-1)'/G;
    NM = EPS*G;
    DM = EPS*NK*G;
    RM = NM-(1/2)*(1-exp(-2*NM));
    RD = DM-(1/2)*(1-exp(-2*DM));

    figure;hold
    plot(NK,sqrt(R)/Ree,'o',NK,sqrt(RD/RM),'-')
    set(gca,'xscale','log');set(gca,'yscale','log')
    xlabel('J/G');ylabel('R_J/R_M')
end

%% Collective conformation tests
if (TEST3)
    % total and partial density
    [PHIA,PHIB] = r_to_phi(r,boxl,DEL,V);
    figure;hist(PHIA+PHIB,50);
    xlabel('\phi_A+phi_B');ylabel('P(\phi)')
end

if (TEST4)
    % radial distrbution function
    [R,g] = gcalc(r,boxl);
    figure;
    plot(R,g(:,1),R,g(:,2),R,g(:,3),R,g(:,4));
    xlabel('r');ylabel('g(r)')
end

if (TEST5)
    % structure factor
    [k,s]=scalc(r,boxl,lksample);
    
    % load analytical theory
    filename = sprintf('sdata/S_EPS%.3fLAM%.2fFA%.2f',EPS,LAM,FA);
    S = load(filename);
    KS = S(:,1);
    SS = 1./(-2*CHI+1./S(:,2));
    
    figure;hold
    plot(k*Ree,s,KS,SS);
    xlabel('q');ylabel('S(q)')
    set(gca,'xscale','log');set(gca,'yscale','log')
end

if (TEST6)
    % chemical potential (widom insertion)
    mu=widom(r,LAM,FA,N,G,CHI,KAP,boxl,DEL,V);
end

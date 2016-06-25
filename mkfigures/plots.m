% plot structure factors
clear; close all
addpath('misc/');
addpath('../utility/');

% simulation parameters
boxl = 20;
Ree = 2;
EPS = 0.01;
LAM = -0.75;
G = 5;

% plot parameters
lksample = 20;
%NREP = [1:2:10,19,29,39];  % number of replicas
%NSNAP = 10:1:10;  % snapshots to average
NREP =40;
NSNAP = 50; 

figure;hold
% plot simulation results
for REP=NREP
    col = (REP-1)/(max(NREP)-1)
    savg = [];
    for SNAP=NSNAP
        r=dlmread(sprintf('../data/r%dv%d',SNAP,REP));
        [k,s]=scalc(r,boxl,lksample);
        savg = [savg,s];
    end
    K = k*Ree; S = mean(savg,2);
    plot(K,S,'color',[col 0 1-col]);
end
[KS,SINV,D2S,ERR]=calcserr(K,S,G)

% plot mean-field theory
CHI = load('../data/cof');
filename = sprintf('../utility/sdata/Seps%.3flam%.2f',EPS,LAM);
spinodal = load('../utility/sdata/chivals');
CHIS = spinodal(spinodal(:,1)==EPS & spinodal(:,2)==LAM,3);
S = load(filename);
KS = S(:,1);
for REP=NREP
    col = (REP-1)/(max(NREP)-1);
    
    if CHI(REP) < CHIS/G
        SS = 1./(-2*CHI(REP)+1./S(:,2));
        plot(KS,SS,'--','color',[col 0 1-col]);
    end
end

% plot processing
set(gca,'xscale','log');set(gca,'yscale','log')

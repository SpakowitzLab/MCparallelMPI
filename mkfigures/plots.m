clear; close all

% plot structure factors
addpath('misc/');
addpath('../utility/');

% simulation parameters
dir = '../../data/';  % data directory
boxl = 20;
Ree = 2;
EPS = 0.01;
LAM = -0.75;
G = 5;

% plot parameters
lksample = 20;
NREP = [1:39];  % number of replicas
NSNAP = 10:1:15;  % snapshots to average

% data to save
KS = zeros(length(NREP),1);
SINV = zeros(length(NREP),1);
D2S = zeros(length(NREP),1);
ERR = zeros(length(NREP),3);

%% Figure 1
figure;hold;cnt=1;
% plot simulation results
for REP=NREP
    col = (REP-1)/(max(NREP)-1)
    savg = [];
    for SNAP=NSNAP
        r=dlmread(strcat(dir,sprintf('r%dv%d',SNAP,REP)));
        [k,s]=scalc(r,boxl,lksample);
        savg = [savg,s];
    end
    K = k*Ree; S = mean(savg,2);
    plot(K,S,'color',[col 0 1-col]);
    
    % fit to Lorentzian
    [KS(cnt),SINV(cnt),D2S(cnt),ERR(cnt,1:3)]=calcserr(K,S,G);
    cnt = cnt+1;
end

% plot mean-field theory
CHI = load(strcat(dir,'cof'));
filename = sprintf('../utility/sdata/Seps%.3flam%.2f',EPS,LAM);
spinodal = load('../utility/sdata/chivals');
CHIS = spinodal(spinodal(:,1)==EPS & spinodal(:,2)==LAM,3);
S = load(filename);
KT = S(:,1);
for REP=NREP
    col = (REP-1)/(max(NREP)-1);
    
    if CHI(REP) < CHIS/G
        ST = 1./(-2*CHI(REP)+1./S(:,2));
        plot(KT,ST,'--','color',[col 0 1-col]);
    end
end

% plot processing
set(gca,'xscale','log');set(gca,'yscale','log')

%% Figure 2
figure;hold;
errorbar(CHI(NREP)*G,SINV,ERR(:,2),'o-')
plot(CHI*G,2*(CHIS/G-CHI))
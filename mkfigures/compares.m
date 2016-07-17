clear;close all

% OPTIONS
PLOTS = 1;      % plot structure factors
PLOTSINV = 1;   % plot inverse peak intensities

% plot parameters
SCALE = 46;     % scaling of simulation structure factor

% add/define paths
dir = '../data/';       % data directory
savedir = 'savedata/';   % save directory

%% Figure 1: structure factors
% load Simulation CHI parameters
G = 5;LAM = 0.0;EPS = 0.04;FA=0.16;
CHI = load(strcat(dir,'chi'));  % adjusted CHI values
CHI = CHI(end,2:end);
PLOTRANGE = unique(round(logspace(0,log10(126),50)));

% load data
rm=32.3348;
sexp = load('savedata/SEXP_PEG30MW1500');
sexp_full = load('savedata/PEG30MW1500.csv');

% load simulation data
figure;hold
NREP = 1:5:30;cnt = 1;
for REP = NREP
    col = (cnt-1)/(length(NREP)-1)
    SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
    ssim = load(strcat(savedir,SAVEFILENAME));
    plot(ssim(PLOTRANGE,1),ssim(PLOTRANGE,2)*SCALE,'o-','linewidth',1.5,'color',[col 0 1-col])
    cnt=cnt+1;
end

NEXP = 3:2:9;cnt = 1;
for ii = NEXP
    col = (cnt-1)/(length(NEXP)-1);
    plot(sexp_full(:,1)*rm,sexp_full(:,ii+1),'linewidth',3,'color',[col 0 1-col])
    cnt=cnt+1;
end

% plot process
set(gca,'xscale','log');set(gca,'yscale','log')

%% Figure 2 (Inverse Peak Intensities)
filename = sprintf('../utility/sdata/S_EPS%.3fLAM%.2fFA%.2f',EPS,LAM,FA);
MF = load(filename);
CHIS = MF(1,1);   %spinodal

figure;hold;set(gca,'fontsize',20)
plot(CHI*G,2*(CHIS/G-CHI),'k-.','linewidth',2)

NREP = 1:30;
for REP = NREP
    SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
    ssim = load(strcat(savedir,SAVEFILENAME));
    plot(CHI(REP)*G,1/max(ssim(:,2)),'o-','linewidth',1.5,'color',[col 0 1-col])
    cnt=cnt+1;
end
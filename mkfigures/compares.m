clear;close all

% OPTIONS
PLOTSIM = 1;      % plot simulation structure factors
PLOTMF = 1;       % plot mean-field structure factors
PLOTEXP = 0;      % plot experiment structure factors
PLOTQS = 1;      % plot inverse peak intensities and peak locations

% plot parameters
NREP1 = [1:9:79];  % simulation plot range for structure factors
NREP2 = [1:9:79];  % simulation plot range for peaks in structure factors
ZEROPK = 0;     % if use zero q as peak location in plot 2
NEXP = 3:2:9;   % experiment plot range

% simulation parameters
EPS = 0.01;  % inter-bead segment rigidity (in unit of 2lp)
LAM = 0.;     % degree of chemical correlation

% simulation constants
G = 5;       % number of beads per monomer
FA = 0.5;   % chemical fraction of A species

% plot options
SCALE = 1/46;     % scaling of simulation structure factor
DISCR = 0;        % if only plot s(k) at selective k values

% add/define paths
%loaddir = '../data/';    % data directory
[pathstr,name,ext] = fileparts(pwd);
title = 'randcopoly';ind = findstr(pathstr, title);
loaddir = strcat('../../../sim-',pathstr(ind:end),'/data/');
savedir = 'savedata/';   % save data directory
figdir = 'figures/';     % save figure directory

%% Figure 1: structure factors
CHI = load(strcat(loaddir,'chi'));  % adjusted CHI values
CHI = CHI(end,2:end);

% load simulation data
figure('Position', [100, 100, 1200, 900]);
hold;set(gca,'fontsize',20);leg = {};

if (PLOTSIM)
    for REP = NREP1
      if max(NREP1) == 1
        col = 1;
      else
        col = (REP-1)/(max(NREP1)-1);
      end
      SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
      ssim = load(strcat(savedir,SAVEFILENAME));

      if (DISCR)
        PLOTRANGE = unique(round(logspace(0,log10(length(ssim)),50)));
        plot(ssim(PLOTRANGE,1),ssim(PLOTRANGE,2),'.-','MarkerSize',15,'linewidth',1.5,'color',[col 0 1-col])
      else
        plot(ssim(:,1),ssim(:,2),'.-','MarkerSize',15,'linewidth',1.5,'color',[col 0 1-col])
      end
      leg{REP} = strcat('\chiG = ', sprintf('%.2f', CHI(REP)*G));
    end
end

if (PLOTMF)
    % Plot mean-field theory
    % load MF results
    filename = sprintf('../utility/sdata/S_EPS%.3fLAM%.2fFA%.2f',EPS,LAM,FA);
    MF = load(filename);
    CHIS = MF(1,1);   %spinodal
    for REP=NREP1
        if max(NREP1)==1
            col = 1;
        else
            col = (REP-1)/(max(NREP1)-1);
        end

        if CHI(REP) < CHIS/G
            K = MF(2:end,1); %wavevectors
            S = 1./(-2*CHI(REP)+EPS*MF(2:end,2));
            plot(K,S,'--','color',[col 0 1-col],'linewidth',2);
        end
    end
end

if (PLOTEXP)
    % load data
    rm=32.3348;
    sexp = load('savedata/SEXP_PEG30MW1500');
    sexp_full = load('savedata/PEG30MW1500.csv');

    cnt = 1;
    for ii = NEXP
      if length(NREP1) == 1
        col = 1;
      else
        col = (cnt-1)/(length(NREP1)-1);
      end
      plot(sexp_full(:,1)*rm,sexp_full(:,ii+1)*SCALE,'linewidth',3,'color',[col 0 1-col])
      cnt=cnt+1;
    end
end

% plot process
legend(leg{NREP1});
axis([0.5,15,1e-2,5e3])
xlabel('R_Mq');ylabel('S(q)');box on
set(gca,'xscale','log');set(gca,'yscale','log')
saveas(gcf,strcat(figdir,'sk'),'epsc')

%% Figure 2 (Inverse Peak Intensities)
if (PLOTQS)
  figure('Position', [100, 500, 1200, 900]);

  subplot(2,1,1);hold;set(gca,'fontsize',20)
  filename = sprintf('../utility/sdata/S_EPS%.3fLAM%.2fFA%.2f',EPS,LAM,FA);
  MF = load(filename);
  CHIS = MF(1,1);   % spinodal
  KS = MF(1,2);     % critical wavemode
  plot(CHI*G,2*(CHIS/G-CHI),'k-.','linewidth',2)

  % find peak location
  QS = zeros(max(NREP2),1);
  for REP = NREP2
    SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
    if (ZEROPK)
      QS(REP) = 1;
    else
      ssim = load(strcat(savedir,SAVEFILENAME));
      QS(REP) = find(ssim(:,2) == max(ssim(:,2)));
    end
  end

  for REP = NREP2
    col = (REP-1)/(max(NREP2)-1);
    SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
    ssim = load(strcat(savedir,SAVEFILENAME));
    plot(CHI(REP)*G,1/ssim(QS(REP),2),'o','linewidth',1.5,'color',[col 0 1-col])
  end
  xlabel('\chiG');ylabel('1/S(q^*)')
  ylim([0,2*CHIS/G]);box on

  subplot(2,1,2);hold;set(gca,'fontsize',20)
  plot(CHI*G,repmat(KS,length(CHI),1),'k-.','linewidth',2)
  for REP = NREP2
    col = (REP-1)/(max(NREP2)-1);
    SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
    ssim = load(strcat(savedir,SAVEFILENAME));
    plot(CHI(REP)*G,ssim(QS(REP),1),'o','linewidth',1.5,'color',[col 0 1-col])
  end
  xlabel('\chiG');ylabel('R_Mq^*');box on
end
saveas(gcf,strcat(figdir,'sinv'),'epsc')

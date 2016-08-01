clear;close all

% OPTIONS
PLOTSIM = 1;    % plot simulation structure factor
PLOTMF = 1;     % plot mean-field structure factor
SAVESIM = 1;    % save simulation structure factor to file

% plot parameters
NREP = [1:60];  % number of replicas
NSNAP = 25:63;  % snapshots to average

% add/define paths
addpath('misc/');
addpath('../utility/');
dir = '../data/';       % data directory
savedir = 'savedata/';   % save directory

% simulation parameters
boxl = 20;   % edge size of simulation
Ree = 2.0;   % average end-to-end distance of a monomer
EPS = 0.01;  % inter-bead segment rigidity (in unit of 2lp)
LAM = 0;     % degree of chemical correlation
G = 5;       % number of beads per monomer
FA = 0.16;   % chemical fraction of A species
lksample = 10;

% load CHI parameters
CHI = load(strcat(dir,'chi'));  % adjusted CHI values
CHI = CHI(end,2:end);

%% (Structure factors)
if (PLOTSIM || PLOTMF)
  figure;hold;set(gca,'fontsize',20);cnt=1;
end

if (PLOTSIM)
    % plot simulation results
    for REP=NREP
      if length(NREP) == 1
        col = 1;
      else
        col = (cnt-1)/(length(NREP)-1);
      end
      savg = [];
      for SNAP=NSNAP
          fprintf('REP = %d, SNAP = %d\n',REP,SNAP)
          r=dlmread(strcat(dir,sprintf('r%dv%d',SNAP,REP)));
          [k,s]=scalc(r,boxl,lksample);
          savg = [savg,s];
      end
      K = k*Ree; S = mean(savg,2);
      plot(K,S,'o-','linewidth',1.5,'color',[col 0 1-col]);

      % save to file
      if (SAVESIM)
          SAVEFILENAME = sprintf('SSIM_CHIG%.3fLAM%.2fEPS%.2fFA%.2f',CHI(REP)*G,LAM,EPS,FA);
          dlmwrite(strcat(savedir,SAVEFILENAME),[K,S]);
      end
        
      % fit to Lorentzian
      %[KS,SINV,D2S,ERR]=calcserr(K,S,G);
        
      % find peak location
      %[pks,locs] = findpeaks(S);
      cnt = cnt+1;
    end
end

if (PLOTMF)
    % Plot mean-field theory
    % load MF results
    filename = sprintf('../utility/sdata/S_EPS%.3fLAM%.2fFA%.2f',EPS,LAM,FA);
    MF = load(filename);
    CHIS = MF(1,1);   %spinodal
    for REP=NREP
        if max(NREP)==1
            col = 1;
        else
            col = (REP-1)/(max(NREP)-1);
        end

        if CHI(REP) < CHIS/G
            K = MF(2:end,1); %wavevectors
            S = 1./(-2*CHI(REP)+EPS*MF(2:end,2));
            plot(K,S,'--','color',[col 0 1-col],'linewidth',2);
        end
    end
end

% plot processing
xlabel('R_Mq');ylabel('S(q)');box on
set(gca,'xscale','log');set(gca,'yscale','log')

path_with=sprintf('ReptAndSwap7_19_24/');

repMIN=1;
repSkip=1;
repMAX=47;
snapMin=26;
snapSkip=1;
snapMax=100;

tset_with=snapMin:snapSkip:snapMax;
repset_with=repMIN:repSkip:repMAX;
nrep_with=length(repset_with);
nt_with=length(tset_with);
Chi_with=zeros(1,nrep_with);
Echi_with=zeros(nt_with,nrep_with);

for rep=repset_with
    out=dlmread(strcat(path_with,sprintf('/out1v%d',rep)),'',1,0);
    Chi_with(1,rep)=out(end,13);
    Echi_with(:,rep)=out(tset_with,8);
end

figure(1); hold on
for tt=1:nt_with
    plot(Chi_with,Echi_with(tt,:)./Chi_with,'k.','MarkerSize',1)
end

x_bar=mean(Echi_with,1)./Chi_with;

plot(Chi_with,x_bar,'r-')
xlabel('\chi')
ylabel('<x_\chi>')

%% Thermodynamic Intigration

F=zeros(nrep_with,1);
for rep=2:nrep_with
    F(rep)=F(rep-1)+0.5*(x_bar(rep)+x_bar(rep-1))*...
                        (Chi_with(rep)-Chi_with(rep-1));
end

figure(2)
plot(Chi_with,F,Chi_with,mean(Echi_with,1))
legend('F','H')
title('Thermodynamic Intigration')
xlabel('\chi')
ylabel('F')

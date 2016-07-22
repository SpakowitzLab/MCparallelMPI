path=sprintf('../data/');
NPT=100;
tstart=1000; % exchange at which to start taking data

data=dlmread(strcat(path,sprintf('chi')));
chi_history=data(:,2:end); % h_A(time,rep)
IND=data(:,1);
N_exchange=length(IND);
nRep=size(chi_history,2);

figure(1)
plot(NPT*(1:N_exchange),chi_history)
xlabel('number of mc moves')
ylabel('h of replica')

chi=chi_history(end,:)';
data=dlmread(strcat(path,sprintf('nodeNumber')));
nodeNumber =data(:,2:end);


nodePath=zeros(size(nodeNumber))*NaN;
for tt=1:N_exchange
    for rep=1:nRep
        nodePath(tt,nodeNumber(tt,rep))=chi_history(tt,rep);
    end
end
size(nodePath)
figure(4)
tskip=25;
repSkip=2;
plot(NPT*(1:tskip:N_exchange),nodePath(1:tskip:end,1:repSkip:end))
xlabel('number of mc moves')
ylabel('h of thread')

figure(3)
plot(nodePath(tstart,:),nodePath(tstart:end,:),'b.')
xlabel('\chi at time 100,000')
ylabel('\chi at later times')

% FieldEnergy=zeros(nRep,1)*NaN;
% for rep=1:nRep
%     data=dlmread(strcat(path,sprintf('out1v%d',rep)),'',1,0);
%     FieldEnergy(rep)=data(end,9);
% end
% figure(5)
% plot(h_A,FieldEnergy,'o')
% xlabel('h')
% ylabel('Field Energy')
% 
% figure(5)
% plot(h_A,FieldEnergy./h_A,'o')
% xlabel('h')
% ylabel('x feild')





function timePlot

%dir='../firstInrData5_30_16/';
dir='../data/';
mu_data=dlmread(strcat(dir,'mu'));
nrep=length(mu_data);

repSkip=5;

leglist=[];

for rep=1:repSkip:nrep
    string=strcat(dir,sprintf('out1v%d',rep));
    data=dlmread(string);
    npts=size(data,1);
    % IND, EELAS(1), EELAS(2), EELAS(3), Eint, EKap,...ECHI, EBind, M
    
    leglist{end+1}=num2str(mu_data(rep));

    figure(1)
    subplot(2,3,1)
    hold on
    plot(data(:,1),data(:,2))
    ylabel('bending Energy')
    
    subplot(2,3,2)
    hold on
    plot(data(:,1),data(:,3))
    ylabel('parallel comp Energy')   
    
    subplot(2,3,3)
    hold on
    plot(data(:,1),data(:,4))
    ylabel('Shear Energy') 
    
    subplot(2,3,4)
    hold on
    plot(data(:,1),data(:,5))
    ylabel('Eint') 
    
%     subplot(2,4,5)
%     hold on
%     plot(data(:,1),data(:,6))
%     ylabel('EKap')     
%     
%     subplot(2,4,6)
%     hold on
%     plot(data(:,1),data(:,7))
%     ylabel('ECHI')   
    
    subplot(2,3,5)
    hold on
    plot(data(:,1),data(:,8))
    ylabel('EBind')   
    
    subplot(2,3,6)
    hold on
    plot(data(:,1),data(:,9))
    ylabel('M')  
    
end
legend(leglist)



end
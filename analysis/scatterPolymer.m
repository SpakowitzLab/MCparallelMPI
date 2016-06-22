
%{
for count=1:10:290
    fName=sprintf('data/r%d',count);
    data=dlmread(fName);
    [rows, cols]=size(data);
    
    mean(data(:,4))
end
%}

    disp('reading params')
    params=dlmread('data/params'); % parameters
    NT=params(1);
    N=params(2);
    NB=params(3);
    NP=params(4);
    disp('NP')
    disp(NP)
    NT=params(5);
    G=params(6);

    L0=params(7)
    CHI0=params(8);
    Fpoly=params(9);
    LBox=params(10);
    disp('LBox')
    disp(LBox)
    EU = params(11);
    EM = params(12)
    HP1_Bind = params(13)
    LK = params(14)
    FRMMETH = params(15);
    F_METH = params(16);
    LAM_METH = params(17);
    
    


% make plot of bound fraction vs. time
boundFraction=0;
if boundFraction
    figure(27)
    disp('Reading methalation')
    u0=dlmread('data/u0'); % u0
    METH=logical(u0(:,5));
    nmeth=sum(METH);
    numeth=length(METH)-nmeth;
    
    f_Meth=zeros(1000,1);
    f_UMeth=zeros(1000,1);
    
    for I=1:5:930
        string=sprintf('data/r%d',I)
        data=dlmread(string);
        msum=sum(METH.*data(:,4));
        wsum=sum((~METH).*data(:,4));
        f_Meth(I)=msum/nmeth;
        f_UMeth(I)=wsum/numeth;
        
    end
    
    plot([1:1000],f_Meth,'.',[1:1000],f_UMeth,'.')
    legend('methalted','unmethalated')
    title(sprintf(...
        'Fractions bound lambda=%5.2f, EU=%5.2f, EM=%5.2f\n HP1 Bind=%5.2f, LK=%5.2f, Fpoly=%5.2f'...
        ,LAM_METH,EU,EM,HP1_Bind,LK,Fpoly))
    xlabel('Save point')
    ylabel('fraction bound')
end

% make plot of "time" history
timeHistory=0;
if (timeHistory)
    disp('Reading time history data')
    out3=dlmread('data/out3');
    figure(8362)
    for j=0:3
        hold all
        subplot(2,2,j+1)
        plot(out3(:,1),out3(:,2+3*j),'o','DisplayName',['window ' num2str(j+1)])
        hold on
        subplot(2,2,j+1)
        plot(out3(:,1),out3(:,3+3*j),'-','DisplayName',['MCAMP ' num2str(j+1)])
        hold on
        subplot(2,2,j+1)
        plot(out3(:,1),out3(:,4+3*j),'-.','DisplayName',['PHIT'  num2str(j+1)])
        hold on
    end
    for j=0:3
        subplot(2,2,j+1)
        legend(gca,'show')
    end
    
    figure(37543)
    out1=dlmread('data/out1');
    hold on
    plot(out1(:,1),out1(:,2),'DisplayName','Bending energy')
    plot(out1(:,1),out1(:,3),'DisplayName','Parallel compression energy')
    plot(out1(:,1),out1(:,4),'DisplayName','Perpendicular compression energy')
    legend(gca,'show')
end


% make plot of particle density w.r.t. position rectangular coordinates
boundryStats=0;
if boundryStats
    spotCheckX=zeros(1001,1);
    disp('Reading methalation')
    u0=dlmread('data/u500'); % u0
    METH=logical(u0(:,5));

    nbins=52;
    bintops=(1:nbins)*LBox;
    binA=zeros(nbins,1);
    binMeth=zeros(nbins,1);
    binTotal=zeros(nbins,1);

    for I=800:5:850
        string=sprintf('data/r%d',I)
        data=dlmread(string);
%        data(:,1:2)=mod(data(:,1:2),LBox); % periodic boundary in x and y

        centers = (LBox/(2*nbins):LBox/nbins:LBox);

        deltaBinTotal = hist(data(:,3),centers);
        binTotal=binTotal+deltaBinTotal';
        deltaBinMeth = hist(data(METH,3),centers);
        binMeth=binMeth+deltaBinMeth';
        deltBinA = hist(data( logical(data(:,4)) ,3),centers);
        binA=binA+deltBinA';

        spotCheckX(I)=data(137,1,1);

    end


    % disp('binTotal')
    % disp(binTotal)
    % disp('binA')
    % disp(binA)
    % disp('binMeth')
    % disp(binMeth)

    figure(1)
    plot(centers,binTotal,centers,binA,centers,binMeth)
    xlabel('z position')
    legend('binTotal','binA','binMeth')

    figure(2)
    plot(centers,binA./binTotal,centers,binMeth./binTotal)
    xlabel('z position')
    legend('fraction of beads A','fraction beads Meth')
    
    figure(4)
    plot(1:1001,spotCheckX,'.')
end

% Make a picture of a sigle snap
MakePicture=0;
if MakePicture

    figure(5)

    disp('csv to matlab...')
    IOstr='data/r11';
    data=dlmread(IOstr);

    [rows, cols]=size(data);

    disp('Applying periodic condition')
    data(:,1:3)=mod(data(:,1:3),LBox);


    disp('coloring...')
    c=zeros(rows,3);
    s=ones(rows,1)*30;
    for row=1:rows
        if data(row,4) % one is read blue is zero
            c(row,:)=[1,0,0]; %red
        else
            c(row,:)=[0,0,1]; %blue
        end
    end

    disp('plotting ...')
    scatter3(data(:,1),data(:,2),data(:,3),s,c,'filled')
    title(IOstr)
    
    xlim([0 LBox])
    zlim([0 LBox])
    ylim([0 LBox])

end

% Plot density in spherical coordinates
PlotRadialDensity=0;
if PlotRadialDensity

    disp('csv to matlab...')
    data=dlmread('data/r409');

    [rows, cols]=size(data);
    
    disp('Reading methalation')
    u0=dlmread('data/u0'); % u0
    METH=logical(u0(:,5));
    
    nMeth=sum(METH);
    nUnMeth=rows-nMeth;
    dataMeth=zeros(nMeth,3);
    cMeth=zeros(nMeth,3);
    dataUnMeth=zeros(nUnMeth,3);
    cUnMeth=zeros(nUnMeth,3);

%    disp('Applying periodic condition')
%    data(:,1:3)=mod(data(:,1:3),LBox);

    nVec=[0,0,0,0]; % UnMeth-Unbound, UnMeth-bound, Meth-Unbound, Meth-Bound 
    disp('coloring...')
    c=zeros(rows,3);
    sMeth=ones(nMeth,1)*30;
    sUnMeth=ones(nUnMeth,1)*30;
    iMeth=1; iUnMeth=1;
    for row=1:rows
        if METH(row)
            dataMeth(iMeth,:)=data(row,1:3);
            if data(row,4) % one is read blue is zero
                cMeth(iMeth,:)=[1,0,0]; %red
                nVec(4)=nVec(4)+1;
            else
                cMeth(iMeth,:)=[0,0,1]; %blue
                nVec(3)=nVec(3)+1;
            end 
            iMeth=iMeth+1;
        else
            dataUnMeth(iUnMeth,:)=data(row,1:3);
            if data(row,4) % one is read blue is zero
                cUnMeth(iUnMeth,:)=[1,0,0]; %red
                nVec(2)=nVec(2)+1;
            else
                cUnMeth(iUnMeth,:)=[0,0,1]; %blue
                nVec(1)=nVec(1)+1;
            end
            iUnMeth=iUnMeth+1;
        end
    end
    
    figure(6)
    disp('plotting ...')
    hold on
    scatter3(dataMeth(:,1),dataMeth(:,2),dataMeth(:,3),sMeth,cMeth,'filled')
    hold all
    scatter3(dataUnMeth(:,1),dataUnMeth(:,2),dataUnMeth(:,3),sUnMeth,cUnMeth,'o')

    xlim([0 LBox])
    zlim([0 LBox])
    ylim([0 LBox])
    
    
    disp('nVec')
    disp(nVec)
    data1=zeros(nVec(1),3);
    data2=zeros(nVec(2),3);
    data3=zeros(nVec(3),3);
    data4=zeros(nVec(4),3);
    ivec=[1,1,1,1];
    for row=1:rows
        if METH(row)
            if data(row,4) % Meth-Bound
                data4(ivec(4),:)=data(row,1:3);
                ivec(4)=ivec(4)+1;
            else % Meth-UnBound
                data3(ivec(3),:)=data(row,1:3);
                ivec(3)=ivec(3)+1;
            end 
        else
            if data(row,4) % UnMeth-Bound
                data2(ivec(2),:)=data(row,1:3);
                ivec(2)=ivec(2)+1;
            else %UnMeth-Unbound
                data1(ivec(1),:)=data(row,1:3);
                ivec(1)=ivec(1)+1;
            end
        end
    end
           
    figure(7)
    R=LBox/2;
    
    nbins=50;
    edges=[0:(nbins+1)].*1/(nbins+1);
    fracV=(((data1(:,1)-R).^2 + (data1(:,2)-R).^2  + (data1(:,3)-R).^2).^(3/2))./R^3;
    h=histogram(fracV,edges);
    density1=h.Values;

    fracV=(((data2(:,1)-R).^2 + (data2(:,2)-R).^2  + (data2(:,3)-R).^2).^(3/2))./R^3;
    h=histogram(fracV,edges);
    density2=h.Values;
    
    fracV=(((data3(:,1)-R).^2 + (data3(:,2)-R).^2  + (data3(:,3)-R).^2).^(3/2))./R^3;
    h=histogram(fracV,edges);
    density3=h.Values;
    
    fracV=(((data4(:,1)-R).^2 + (data4(:,2)-R).^2  + (data4(:,3)-R).^2).^(3/2))./R^3;
    h=histogram(fracV,edges);
    density4=h.Values;
    xpts=0.5*(edges(1:end-1)+edges(2:end));
    
    densityTotal=density1+density2+density3+density4;
    
    figure(8)
    hold all
    plot(xpts,density1,'.');
    plot(xpts,density2,'o');
    plot(xpts,density3,'+');
    plot(xpts,density4,'*');
    legend('UnMeth-Unbound','UnMeth-Bound','Meth-Unbound','Meth-Bound')
    xlabel('r^3 / R^3')
    ylabel('density')
    set(gca,'YTick',[])
    
    figure(9)
    plot(xpts,(density1+density3)./(density2+density4),'o')
    title('denisty unbound/density bound')
    xlabel('r^3 / R^3')
    ylabel('fraction ')
    
    
    
    figure(10)
    plot(xpts.^0.33333,densityTotal*nbins/(8000*4*pi()/3),'o')
    hold on
    Rvec=(0:200)/200;
    Yvec=(sin(pi()*Rvec).^2)./(2*pi()*Rvec.^2);
    plot(Rvec,Yvec,'+')
    
    xlabel('r/R')
    ylabel('relitiveDensity')
    xlim([0,1])
  
    disp('total number')
    disp(sum(densityTotal))

end

makeRandom=0;
if makeRandom
    npts=NT;
    data=zeros(npts,3);
    R=LBox/2;
    for ii=1:npts
        stop=1;
        while(stop)
            t=rand([1,3])*LBox;
            if ((norm(t-[R,R,R]))<R)
                data(ii,:)=t;
                stop=0;
            end
        end
    end
    fprintf('R=%f, LBOX=%f\n',R,LBox)
end

particleTracks=0;
if particleTracks
    max=200;
    step=5;
    npts=max/step +1;
    history=zeros(1,npts);
    n=0;
    for I=0:step:max
        string=sprintf('data/r%d',I);
        data=dlmread(string);
        n=n+1;
    end    
end

freeDensityCompair=1;
if freeDensityCompair
    nbins=40;
    overset=10
    density=zeros(1,nbins+overset);
    n=0;
    R=LBox/2;
    edges=(0:nbins+overset).*1/(nbins);
    xpts=0.5*(edges(1:end-1)+edges(2:end));
    disp('loading for freeDensityCompair')
    for I=450:5:499
        string=sprintf('data/r%d',I);
        data=dlmread(string);
        fracV=(((data(:,1)-R).^2 + (data(:,2)-R).^2  + (data(:,3)-R).^2).^(3/2))./R^3;
        h=histogram(fracV,edges);
        density=density+h.Values;
        n=n+1;
    end
    density=density/n;
    figure(121)
    plot(xpts.^0.33333,density*nbins/(NT*4*pi()/3),'o')
    hold on
    Rvec=(0:200)/200;
    Yvec=(sin(pi()*Rvec).^2)./(2*pi()*Rvec.^2);
    plot(Rvec,Yvec,'-k')
    title(sprintf('L0=%f, LK=%f, R=%f, NT=%d',L0,LK,R,NT));
    xlabel('r/R')
    ylabel('relative density')
    %disp('Density')
    %disp(density)
    disp('Average x,y,z')
    disp([mean(data(:,1)),mean(data(:,2)),mean(data(:,3))])
    
end

sloshPlot=1
if sloshPlot
    nbins=20;
    density=zeros(1,nbins);
    R=LBox/2;
    edges=(0:nbins).*1/(nbins);
    xpts=0.5*(edges(1:end-1)+edges(2:end));
    disp('loading for slosh plot')
    colors=[1 0 0;
            0.9 0 0.1;
            0.8 0 0.2;
            0.7 0 0.3;
            0.6 0 0.4;
            0.5 0 0.5;
            0.4 0 0.6;
            0.3 0 0.7;
            0.2 0 0.8;
            0.1 0 0.9;
            0 0 1];
    n=1;
    for I=[0,2,4,8,16,32,64,128,256]
        %string=sprintf('slochData/r%d',I);
        string=sprintf('data/r%d',I);
        data=dlmread(string);
        fracV=(((data(:,1)-R).^2 + (data(:,2)-R).^2  + (data(:,3)-R).^2).^(3/2))./R^3;
        figure(128)
        h=histogram(fracV,edges);
        density=h.Values;
        figure(21)
        hold on
        plot(xpts.^0.33333,density*nbins/(NT*4*pi()/3),'color',colors(n,:))
        n=n+1;
    end
    figure(21)
    hold on
    Rvec=(0:200)/200;
    Yvec=(sin(pi()*Rvec).^2)./(2*pi()*Rvec.^2);
    plot(Rvec,Yvec,'-k')
    title(sprintf('L0=%f, LK=%f, R=%f, NT=%d',L0,LK,R,NT));
    xlabel('r/R')
    ylabel('relative density')
    legend('0','2','4','8','16','32','64','128','256','GC')
    %disp('Density')
    %disp(density)
    disp('Average x,y,z')
    disp([mean(data(:,1)),mean(data(:,2)),mean(data(:,3))]) 
end


plotSteps=0;
if plotSteps
    string=sprintf('data/r%d',100);
    data=dlmread(string); 
    size(data)
    
    nstep=100;
    steps=data(nstep+1:nstep:end,1:3)-data(1:nstep:end-nstep,1:3);
    drs=(steps(:,1).^2 + steps(:,2).^2 + steps(:,3).^2).^0.5;
    figure(12)
    histogram(drs,100)
    title('distance to bead nstep away')
    fprintf('mean displacement after %d beads is %f\n',nstep,mean(drs))
    fprintf('r giration=mean*(NT/(6*nstep))^0.5=%f\n',mean(drs)*(NT/(6*nstep))^0.5)
end
   
    
    


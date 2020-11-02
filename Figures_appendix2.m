close all
clear all

load('./Data/tables.mat')

eq = readtable('./Data/relacio_aga.xlsx');
prop = xlsread('./Figures&Tables/Table_incidence_T2.xlsx','Numbers','L:L');
prop = [prop(2:end);prop(1)];
prop = 1./prop;

time = datenum('01/06/2020','dd/mm/yyyy'):datenum('03/09/2020','dd/mm/yyyy');
u = unique(T0827.id_AGA);

pop = 0.*u;
for i = 1:length(u)+1
    if i<=length(u)
        pop(i) = unique(T0827.popAGA(T0827.id_AGA==u(i)));
    else
        pop(i) = sum(pop);
    end
end

A14_Pc = zeros(length(time),length(u)+1);
A14_CC = zeros(length(time),length(u)+1);
iA14_Pc = zeros(length(time),length(u)+1);
iA14_CC = zeros(length(time),length(u)+1);

cas_Pc = zeros(length(time),length(u)+1);
cas_CC = zeros(length(time),length(u)+1);
icas_Pc = zeros(length(time),length(u)+1);
icas_CC = zeros(length(time),length(u)+1);

cas7_Pc = zeros(length(time),length(u)+1);
cas7_CC = zeros(length(time),length(u)+1);
icas7_Pc = zeros(length(time),length(u)+1);
icas7_CC = zeros(length(time),length(u)+1);

rho_Pc = zeros(length(time),length(u)+1);
rho_CC = zeros(length(time),length(u)+1);
irho_Pc = zeros(length(time),length(u)+1);
irho_CC = zeros(length(time),length(u)+1);

T = T0913;

for i = 1:length(u)+1
    
    for k = 1:length(time)
        
        if i <= length(u)
            cas_Pc(k,i) = T.PCRcas(T.id_AGA==u(i) & datenum(T.data)==time(k));
            cas_CC(k,i) = T.ClinicsCovid(T.id_AGA==u(i) & datenum(T.data)==time(k));
        else
            cas_Pc(k,i) = sum(T.PCRcas(datenum(T.data)==time(k)));
            cas_CC(k,i) = sum(T.ClinicsCovid(datenum(T.data)==time(k)));
        end
        
        [phat, pci] = binofit(cas_Pc(k,i),pop(i),0.05);
        icas_Pc(k,i) = mean(abs(pci-phat)*pop(i));
        [phat, pci] = binofit(cas_CC(k,i),pop(i),0.05);
        icas_CC(k,i) = mean(abs(pci-phat)*pop(i));
        
        if k>6
            cas7_Pc(k,i) = 1/7*sum(cas_Pc(k-6:k,i));
            cas7_CC(k,i) = 1/7*sum(cas_CC(k-6:k,i));
            icas7_Pc(k,i) = 1/7*sqrt(sum(icas_Pc(k-6:k,i).^2));
            icas7_CC(k,i) = 1/7*sqrt(sum(icas_CC(k-6:k,i).^2));
        end
        
    end
    
    id = 8:length(time);
    rho_Pc(id,i) = (cas7_Pc(id,i)+cas7_Pc(id-1,i)+cas7_Pc(id-2,i))./max(cas7_Pc(id-5,i)+cas7_Pc(id-6,i)+cas7_Pc(id-7,i),1);
    rho_CC(id,i) = (cas7_CC(id,i)+cas7_CC(id-1,i)+cas7_CC(id-2,i))./max(cas7_CC(id-5,i)+cas7_CC(id-6,i)+cas7_CC(id-7,i),1);
    
    for k = 8:length(time)
        k2 = 1/sum(cas7_Pc(k-7:k-5,i));
        irho_Pc(k,i) = k2*sqrt(sum((icas7_Pc(k-2:k,i)).^2)+rho_Pc(k,i)^2*sum((icas7_Pc(k-7:k-5,i)).^2));
        k2 = 1/sum(cas7_CC(k-7:k-5,i));
        irho_CC(k,i) = k2*sqrt(sum((icas7_CC(k-2:k,i)).^2)+rho_CC(k,i)^2*sum((icas7_CC(k-7:k-5,i)).^2));
    end
    
    acc = cumsum(cas7_Pc(:,i));
    a14 = 0.*acc;
    a14(1:14) = acc(1:14);
    a14(15:end) = acc(15:end)-acc(1:end-14);
    A14_Pc(:,i) = a14/pop(i)*1e5;
    acc = cumsum(cas7_CC(:,i));
    a14 = 0.*acc;
    a14(1:14) = acc(1:14);
    a14(15:end) = acc(15:end)-acc(1:end-14);
    A14_CC(:,i) = a14/pop(i)*1e5;
    
    
    for k = 14:length(time)
        iA14_Pc(k,i) = sqrt(sum((icas7_Pc(k-13:k,i)).^2))/pop(i)*1e5;
        iA14_CC(k,i) = sqrt(sum((icas7_CC(k-13:k,i)).^2))/pop(i)*1e5;
    end
    
end

EPG_Pc = A14_Pc.*rho_Pc;
EPG_CC = A14_CC.*rho_CC;

iEPG_Pc = sqrt(A14_Pc.^2.*(irho_Pc.^2) + rho_Pc.^2.*(iA14_Pc.^2));
iEPG_CC = sqrt(A14_CC.^2.*(irho_CC.^2) + rho_CC.^2.*(iA14_CC.^2));

col = lines(4);
u = [u;0];

Name = cell(length(u),1);
for k = 1:length(u)-1
    Name(k) = eq.aga(eq.id_aga==u(k));
end
Name{end} = 'CATALUNYA';

time(1:31)=[];
for k = 1:2
    if k==2
        s1 = 'i';
    else
        s1 = '';
    end
    for i = 1:5
        switch i
            case 1
                s2 = 'cas';
            case 2
                s2 = 'cas7';
            case 3
                s2 = 'A14';
            case 4
                s2 = 'EPG';
            case 5
                s2 = 'rho';
        end
        for j = 1:2
            if j == 1
                s3 = 'Pc';
            else
                s3 = 'CC';
            end
            eval([s1 s2 '_' s3 '(1:31,:) = [];'])
        end
    end
end

%% FIGURA 2A

for ID = 1:43
    
    uy = 0.1;
    sy = 0.04;
    dy = 0.02;
    ux = 0.08;
    sx = 0.08;
    dx = 0.02;
    DY = (1-uy-sy-dy)/2;
    DX = (1-ux-sx-dx)/2;
    
    f2a = figure(1);
    clf
    f2a.Position = [33   371   861   591];
    
    ax = axes('Units','normalized','Position',[ux uy+DY+sy DX DY]);
    
    plot(time,cas7_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,prop(ID)*cas7_CC(:,ID),'linewidth',2,'Color',col(2,:))
    YL = ax.YLim;
    YL(1) = 0;
    
    opt = {'Pc','CC'};
    for s = 1:length(opt)
        eval(['x = [cas7_' opt{s} '(:,ID)+icas7_' opt{s} '(:,ID); flip(cas7_' opt{s} '(:,ID)-icas7_' opt{s} '(:,ID))];'])
        x = x';
        if s==2
            x = x*prop(ID);
        end
        x(x<YL(1)) = YL(1);
        x(x>YL(2)) = YL(2);
        fill([time flip(time)],x,'r','FaceColor',0.5+0.5.*col(s,:),'EdgeColor','none','FaceAlpha',0.5)
    end
    
    plot(time,cas7_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,prop(ID)*cas7_CC(:,ID),'linewidth',2,'Color',col(2,:))
    
    legend('PCR','EMR','Location','best')
    
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    ax.TickDir = 'out';
    
    ax.XLim = [min(time) max(time)];
    ax.YLim = YL;
    xt = min(time):7:max(time);
    ax.XTick = xt;
    XT = cell(length(xt),1);
    for kk = 1:length(xt)
        XT{kk} = datestr(xt(kk),'dd/mm');
    end
    ax.XTickLabel = [];%XT;
    xtickangle(45)
    ylabel('Average cases last 7 days')
    text(0.01,1.02,'A','Units','normalized','FontSize',24,'HorizontalAlignment','left','VerticalAlignment','top')
    
    
    ax = axes('Units','normalized','Position',[ux+DX+sx uy+DY+sy DX DY]);
    
    plot(time,rho_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,rho_CC(:,ID),'linewidth',2,'Color',col(2,:))
    plot(time,0.*time+1,'-k')
    
    opt = {'Pc','CC'};
    for s = 1:length(opt)
        eval(['x = [rho_' opt{s} '(:,ID)+irho_' opt{s} '(:,ID); flip(rho_' opt{s} '(:,ID)-irho_' opt{s} '(:,ID))];'])
        x = x';
        x(x<0) = 0;
        x(x>5) = 5;
        fill([time flip(time)],x,'r','FaceColor',0.5+0.5.*col(s,:),'EdgeColor','none','FaceAlpha',0.5)
    end
    
    plot(time,rho_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,rho_CC(:,ID),'linewidth',2,'Color',col(2,:))
    
    %legend('PCR','EMR','Location','best')
    
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    ax.TickDir = 'out';
    
    axis([min(time) max(time) 0 3])
    xt = min(time):7:max(time);
    ax.XTick = xt;
    XT = cell(length(xt),1);
    for kk = 1:length(xt)
        XT{kk} = datestr(xt(kk),'dd/mm');
    end
    ax.XTickLabel = [];%XT;
    xtickangle(45)
    ylabel('Empiric propagation')
    text(0.01,1.02,'B','Units','normalized','FontSize',24,'HorizontalAlignment','left','VerticalAlignment','top')
    
    ax = axes('Units','normalized','Position',[ux uy DX DY]);
    
    plot(time,A14_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,prop(ID)*A14_CC(:,ID),'linewidth',2,'Color',col(2,:))
    YL = ax.YLim;
    YL(1) = 0;
    
    opt = {'Pc','CC'};
    for s = 1:length(opt)
        eval(['x = [A14_' opt{s} '(:,ID)+iA14_' opt{s} '(:,ID); flip(A14_' opt{s} '(:,ID)-iA14_' opt{s} '(:,ID))];'])
        x = x';
        if s==2
            x = x*prop(ID);
        end
        x(x<YL(1)) = YL(1);
        x(x>YL(2)) = YL(2);
        fill([time flip(time)],x,'r','FaceColor',0.5+0.5.*col(s,:),'EdgeColor','none','FaceAlpha',0.5)
    end
    
    plot(time,A14_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,prop(ID)*A14_CC(:,ID),'linewidth',2,'Color',col(2,:))
    
    %legend('PCR','EMR','Location','best')
    
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    ax.TickDir = 'out';
    
    ax.XLim = [min(time) max(time)];
    ax.YLim = YL;
    xt = min(time):7:max(time);
    ax.XTick = xt;
    XT = cell(length(xt),1);
    for kk = 1:length(xt)
        XT{kk} = datestr(xt(kk),'dd/mm');
    end
    ax.XTickLabel = XT;
    xtickangle(45)
    ylabel('A14')
    text(0.01,1.02,'C','Units','normalized','FontSize',24,'HorizontalAlignment','left','VerticalAlignment','top')
    
    ax = axes('Units','normalized','Position',[ux+DX+sx uy DX DY]);
    
    plot(time,EPG_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,prop(ID)*EPG_CC(:,ID),'linewidth',2,'Color',col(2,:))
    YL = ax.YLim;
    YL(1) = 0;
    
    opt = {'Pc','CC'};
    for s = 1:length(opt)
        eval(['x = [EPG_' opt{s} '(:,ID)+iEPG_' opt{s} '(:,ID); flip(EPG_' opt{s} '(:,ID)-iEPG_' opt{s} '(:,ID))];'])
        x = x';
        if s==2
            x = x*prop(ID);
        end
        x(x<YL(1)) = YL(1);
        x(x>YL(2)) = YL(2);
        fill([time flip(time)],x,'r','FaceColor',0.5+0.5.*col(s,:),'EdgeColor','none','FaceAlpha',0.5)
    end
    
    plot(time,EPG_Pc(:,ID),'linewidth',2,'Color',col(1,:))
    hold on
    plot(time,prop(ID)*EPG_CC(:,ID),'linewidth',2,'Color',col(2,:))
    
    %legend('PCR','EMR','Location','best')
    
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    ax.TickDir = 'out';
    
    ax.XLim = [min(time) max(time)];
    ax.YLim = YL;
    xt = min(time):7:max(time);
    ax.XTick = xt;
    XT = cell(length(xt),1);
    for kk = 1:length(xt)
        XT{kk} = datestr(xt(kk),'dd/mm');
    end
    ax.XTickLabel = XT;
    xtickangle(45)
    ylabel('EPG')
    text(0.01,1.02,'D','Units','normalized','FontSize',24,'HorizontalAlignment','left','VerticalAlignment','top')
    
    print(f2a,['./Figures&Tables/Appendix2/AGA_' num2str(ID,'%02u') '.png'],'-dpng','-r600')
    
end
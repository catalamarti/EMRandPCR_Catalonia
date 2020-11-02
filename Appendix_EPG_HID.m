close all
clear all

load('./Data/tables.mat')
TG = readtable('./Data/dades_risc_gravetat_2020-09-28.txt','Delimiter','@');

eq = readtable('./Data/relacio_aga.xlsx');
prop = xlsread('./Figures&Tables/Table_incidence_T2.xlsx','Numbers','L:L');
prop = [prop(2:end);prop(1)];
prop = 1./prop;

time = datenum('01/03/2020','dd/mm/yyyy'):datenum('03/09/2020','dd/mm/yyyy');
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

Hosp = zeros(length(time),length(u)+1);
iHosp = zeros(length(time),length(u)+1);
UCIs = zeros(length(time),length(u)+1);
iUCIs = zeros(length(time),length(u)+1);
morts = zeros(length(time),length(u)+1);
imorts = zeros(length(time),length(u)+1,2);

Hosp7 = zeros(length(time),length(u)+1);
iHosp7 = zeros(length(time),length(u)+1);
UCIs7 = zeros(length(time),length(u)+1);
iUCIs7 = zeros(length(time),length(u)+1);
morts7 = zeros(length(time),length(u)+1);
imorts7 = zeros(length(time),length(u)+1,2);

T = T0913;

for i = 1:length(u)+1
    
    for k = 1:length(time)
        
        if i <= length(u)
            cas_Pc(k,i) = T.PCRcas(T.id_AGA==u(i) & datenum(T.data)==time(k));
            cas_CC(k,i) = T.ClinicsCovid(T.id_AGA==u(i) & datenum(T.data)==time(k));
            Hosp(k,i) = TG.Var3(TG.Var2==u(i) & datenum(TG.Var1)==time(k));
            UCIs(k,i) = TG.Var4(TG.Var2==u(i) & datenum(TG.Var1)==time(k));
            morts(k,i) = TG.Var5(TG.Var2==u(i) & datenum(TG.Var1)==time(k));
        else
            cas_Pc(k,i) = sum(T.PCRcas(datenum(T.data)==time(k)));
            cas_CC(k,i) = sum(T.ClinicsCovid(datenum(T.data)==time(k)));
            Hosp(k,i) = sum(TG.Var3(datenum(TG.Var1)==time(k)));
            UCIs(k,i) = sum(TG.Var4(datenum(TG.Var1)==time(k)));
            morts(k,i) = sum(TG.Var5(datenum(TG.Var1)==time(k)));
        end
        
        [phat, pci] = binofit(cas_Pc(k,i),pop(i),0.05);
        icas_Pc(k,i) = mean(abs(pci-phat)*pop(i));
        [phat, pci] = binofit(cas_CC(k,i),pop(i),0.05);
        icas_CC(k,i) = mean(abs(pci-phat)*pop(i));
        [phat, pci] = binofit(Hosp(k,i),pop(i),0.05);
        iHosp(k,i) = mean(abs(pci-phat)*pop(i));
        [phat, pci] = binofit(UCIs(k,i),pop(i),0.05);
        iUCIs(k,i) = mean(abs(pci-phat)*pop(i));
        [phat, pci] = binofit(morts(k,i),pop(i),0.05);
        imorts(k,i,1) = (phat-pci(1))*pop(i);
        imorts(k,i,2) = (pci(2)-phat)*pop(i);
        
        Hosp(k,i) = Hosp(k,i)/pop(i)*1e5;
        iHosp(k,i) = iHosp(k,i)/pop(i)*1e5;
        morts(k,i) = morts(k,i)/pop(i)*1e5;
        imorts(k,i,:) = imorts(k,i,:)/pop(i)*1e5;
        UCIs(k,i) = UCIs(k,i)/pop(i)*1e5;
        iUCIs(k,i) = iUCIs(k,i)/pop(i)*1e5;
        
        if k>6
            cas7_Pc(k,i) = 1/7*sum(cas_Pc(k-6:k,i));
            cas7_CC(k,i) = 1/7*sum(cas_CC(k-6:k,i));
            icas7_Pc(k,i) = 1/7*sqrt(sum(icas_Pc(k-6:k,i).^2));
            icas7_CC(k,i) = 1/7*sqrt(sum(icas_CC(k-6:k,i).^2));
            Hosp7(k,i) = 1/7*sum(Hosp(k-6:k,i));
            iHosp7(k,i) = 1/7*sqrt(sum(iHosp(k-6:k,i).^2));
            UCIs7(k,i) = 1/7*sum(UCIs(k-6:k,i));
            iUCIs7(k,i) = 1/7*sqrt(sum(iUCIs(k-6:k,i).^2));
            morts7(k,i) = 1/7*sum(morts(k-6:k,i));
            imorts7(k,i,1) = 1/7*sqrt(sum(imorts(k-6:k,i,1).^2));
            imorts7(k,i,2) = 1/7*sqrt(sum(imorts(k-6:k,i,2).^2));
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

EPG_CC = (prop.*EPG_CC')';
iEPG_CC = (prop.*iEPG_CC')';

u = [u;0];

Name = cell(length(u),1);
for k = 1:length(u)-1
    Name(k) = eq.aga(eq.id_aga==u(k));
end
Name{end} = 'CATALUNYA';

%%

for id = 1:44
    
    f = figure(1);
    clf
    
    xt = datenum('01/07/2020','dd/mm/yyyy'):5:datenum('31/08/2020','dd/mm/yyyy');
    XT = cell(1,length(xt));
    for k = 1:length(xt)
        XT{k} = datestr(xt(k),'dd/mm');
    end
    
    yyaxis left
    plot(time,EPG_Pc(:,id),'-k','linewidth',1.5)
    hold on
    plot(time,EPG_CC(:,id),'-','linewidth',1.5,'Color',[0.6 0.1 0])
    
    ax = gca;
    ax.XLim = [min(xt) max(xt)];
    ax.XTick = xt;
    ax.XTickLabel = XT;
    xtickangle(45)
    ax.TickDir = 'out';
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    ax.YColor = [1 1 1]*0.15;
    ylabel('EPG')
    YL1 = ax.YLim;
    
    nu = 3;
    nm = 30;
    
    yyaxis right
    plot(time,nm*morts7(:,id),'-b','linewidth',1.5)
    hold on
    plot(time,Hosp7(:,id),'-g','linewidth',1.5)
    plot(time,nu*UCIs7(:,id),'-r','linewidth',1.5)
    
    ax = gca;
    ax.XLim = [min(xt) max(xt)+1];
    ax.XTick = xt;
    ax.XTickLabel = XT;
    xtickangle(45)
    ax.TickDir = 'out';
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    ax.YColor = [1 1 1]*0.15;
    ylabel('Deaths-Hospitalizated-ICUs per 10^5 inh.')
    YL2 = ax.YLim;
    m = max(YL1(2)/250,YL2(2)/15);
    ax.YLim = [0 15*m];
    
    yyaxis left
    p1=plot(time,EPG_Pc(:,id),'-k','linewidth',1.5);
    p2=plot(time,EPG_CC(:,id),'-','linewidth',1.5,'Color',[0.6 0.1 0]);
    ax.YLim = [0 250*m];
    
    title(Name{id})
    
    l1 = legend([p1 p2],{'PCR EPG','EMR EPG'},'Location','northwest');
    axes('Units','normalized','Position',[-1 -1 0.1 0.1])
    p3=plot(time,nm*morts7(:,id),'-b','linewidth',1.5);
    hold on
    p4=plot(time,Hosp7(:,id),'-g','linewidth',1.5);
    p5=plot(time,nu*UCIs7(:,id),'-r','linewidth',1.5);
    l2 = legend([p3 p4 p5],{'Deaths (x30)','Hospitalized','ICUs (x3)'});
    l2.Position(1:2) = [0.8857-l2.Position(3) 0.9-l2.Position(4)];
    l2.FontSize = 9;
    l2.LineWidth = 1.5;
    
    print(f,['./Figures&Tables/Appendix_EPG_HID/EPG_HID_' num2str(id,'%02u') '.png'],'-dpng','-r600')
    
end
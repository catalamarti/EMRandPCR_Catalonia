close all
clear all

load('./Data/tables.mat')
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

A14_PCR = zeros(length(time),length(u)+1);
A14_EMR = zeros(length(time),length(u)+1);
cas_PCR = zeros(length(time),length(u)+1);
cas_EMR = zeros(length(time),length(u)+1);
cas7_PCR = zeros(length(time),length(u)+1);
cas7_EMR = zeros(length(time),length(u)+1);
rho_PCR = zeros(length(time),length(u)+1);
rho_EMR = zeros(length(time),length(u)+1);
iA14_PCR = zeros(length(time),length(u)+1);
iA14_EMR = zeros(length(time),length(u)+1);
icas_PCR = zeros(length(time),length(u)+1);
icas_EMR = zeros(length(time),length(u)+1);
icas7_PCR = zeros(length(time),length(u)+1);
icas7_EMR = zeros(length(time),length(u)+1);
irho_PCR = zeros(length(time),length(u)+1);
irho_EMR = zeros(length(time),length(u)+1);

T = T0913;

for i = 1:length(u)+1
    
    for k = 1:length(time)
        
        if i <= length(u)
            cas_PCR(k,i) = T.PCRcas(T.id_AGA==u(i) & datenum(T.data)==time(k));
            cas_EMR(k,i) = T.ClinicsCovid(T.id_AGA==u(i) & datenum(T.data)==time(k));
        else
            cas_PCR(k,i) = sum(T.PCRcas(datenum(T.data)==time(k)));
            cas_EMR(k,i) = sum(T.ClinicsCovid(datenum(T.data)==time(k)));
        end
        icas_PCR(k,i) = 1.96*sqrt(cas_PCR(k,i));
        icas_EMR(k,i) = 1.96*sqrt(cas_EMR(k,i));
        
        if k>6
            cas7_PCR(k,i) = 1/7*sum(cas_PCR(k-6:k,i));
            cas7_EMR(k,i) = 1/7*sum(cas_EMR(k-6:k,i));
            icas7_PCR(k,i) = 1/7*sqrt(sum(icas_PCR(k-6:k,i).^2));
            icas7_EMR(k,i) = 1/7*sqrt(sum(icas_EMR(k-6:k,i).^2));
        end
        
    end
    
    id = 8:length(time);
    rho_PCR(id,i) = (cas7_PCR(id,i)+cas7_PCR(id-1,i)+cas7_PCR(id-2,i))./max(cas7_PCR(id-5,i)+cas7_PCR(id-6,i)+cas7_PCR(id-7,i),1);
    rho_EMR(id,i) = (cas7_EMR(id,i)+cas7_EMR(id-1,i)+cas7_EMR(id-2,i))./max(cas7_EMR(id-5,i)+cas7_EMR(id-6,i)+cas7_EMR(id-7,i),1);
    for k = 8:length(time)
        k2 = 1/sum(cas7_PCR(k-7:k-5,i));
        irho_PCR(k,i) = k2*sqrt(sum((icas7_PCR(k-2:k,i)).^2)+rho_PCR(k,i)^2*sum((icas7_PCR(k-7:k-5,i)).^2));
        k2 = 1/sum(cas7_EMR(k-7:k-5,i));
        irho_EMR(k,i) = k2*sqrt(sum((icas7_EMR(k-2:k,i)).^2)+rho_EMR(k,i)^2*sum((icas7_EMR(k-7:k-5,i)).^2));
    end
    
    acc = cumsum(cas7_PCR(:,i));
    a14 = 0.*acc;
    a14(1:14) = acc(1:14);
    a14(15:end) = acc(15:end)-acc(1:end-14);
    A14_PCR(:,i) = a14/pop(i)*1e5;
    acc = cumsum(cas7_EMR(:,i));
    a14 = 0.*acc;
    a14(1:14) = acc(1:14);
    a14(15:end) = acc(15:end)-acc(1:end-14);
    A14_EMR(:,i) = a14/pop(i)*1e5;
    
    for k = 14:length(time)
        iA14_PCR(k,i) = sqrt(sum((icas7_PCR(k-13:k,i)).^2))/pop(i)*1e5;
        iA14_EMR(k,i) = sqrt(sum((icas7_EMR(k-13:k,i)).^2))/pop(i)*1e5;
    end
    
end

EPG_PCR = A14_PCR.*rho_PCR;
EPG_EMR = A14_EMR.*rho_EMR;

iEPG_PCR = sqrt(A14_PCR.^2.*(irho_PCR.^2) + rho_PCR.^2.*(iA14_PCR.^2));
iEPG_EMR = sqrt(A14_EMR.^2.*(irho_EMR.^2) + rho_EMR.^2.*(iA14_EMR.^2));

EPG_EMR = (prop.*EPG_EMR')';
iEPG_EMR = (prop.*iEPG_EMR')';

EPG_PCR_JA = EPG_PCR(month(time)==7 | month(time)==8,:);
EPG_EMR_JA = EPG_EMR(month(time)==7 | month(time)==8,:);
iEPG_PCR_JA = iEPG_PCR(month(time)==7 | month(time)==8,:);
iEPG_EMR_JA = iEPG_EMR(month(time)==7 | month(time)==8,:);

%%
f = figure(1);
clf

di = abs(EPG_PCR_JA(:)-EPG_EMR_JA(:));
di(di>100) = 100;

scatter(EPG_PCR_JA(:),EPG_EMR_JA(:),10,di,'filled')
axis equal
ax = gca;
ax.XLim = [0 200];
ax.YLim = [0 200];
ax.YTick = ax.XTick;
ax.LineWidth = 1.5;
ax.TickDir='out';
xlabel('EPG PCR')
ylabel('EPG EMR')

%print(f,'./Figures&Tables/Figure_EPG_SFZ.png','-dpng','-r600')

for ID = 1:44
    
    f = figure(1);
    clf
    
    di = abs(EPG_PCR_JA(:,ID)-EPG_EMR_JA(:,ID));
    di(di>100) = 100;
    
    rl_PCR = risklevel(EPG_PCR_JA(:,ID));
    rl_EMR = risklevel(EPG_EMR_JA(:,ID));
    over = overlap(EPG_PCR_JA(:,ID),iEPG_PCR_JA(:,ID),EPG_EMR_JA(:,ID),iEPG_EMR_JA(:,ID));
    c = over | rl_PCR==rl_EMR;
    
    plot([0 200],[0 200],'-k')
    hold on
    fill([0 30 30 0 0],[0 0 30 30 0],'r','FaceColor',[0.4 1 0.4],'EdgeColor','none','FaceAlpha',0.5)
    fill([30 70 70 30 30],[30 30 70 70 30],'r','FaceColor',[1 1 0.4],'EdgeColor','none','FaceAlpha',0.5)
    fill([70 100 100 70 70],[70 70 100 100 70],'r','FaceColor',[1 0.7 0.4],'EdgeColor','none','FaceAlpha',0.5)
    fill([100 200 200 100 100],[100 100 200 200 100],'r','FaceColor',[1 0.4 0.4],'EdgeColor','none','FaceAlpha',0.5)
    plot([0 200],[0 200],'-k')
    errorbar(EPG_PCR_JA(c==0,ID),EPG_EMR_JA(c==0,ID),iEPG_EMR_JA(c==0,ID),iEPG_EMR_JA(c==0,ID),iEPG_PCR_JA(c==0,ID),iEPG_PCR_JA(c==0,ID),'.k','CapSize',0,'Color',col_di(0))
    errorbar(EPG_PCR_JA(c==1,ID),EPG_EMR_JA(c==1,ID),iEPG_EMR_JA(c==1,ID),iEPG_EMR_JA(c==1,ID),iEPG_PCR_JA(c==1,ID),iEPG_PCR_JA(c==1,ID),'.k','CapSize',0,'Color',col_di(1))
    plot([0 200],[0 200],'-k')
    scatter(EPG_PCR_JA(:,ID),EPG_EMR_JA(:,ID),10,color_error(c),'filled')
    axis equal
    ax = gca;
    ax.Box = 'off';
    ax.XLim = [0 200];
    ax.YLim = [0 200];
    ax.YTick = ax.XTick;
    ax.LineWidth = 1.5;
    ax.TickDir='out';
    xlabel('EPG PCR')
    ylabel('EPG EMR')
    print(f,['./Figures&Tables/Appendix_EPG_scatter/Figure_EPG_SFZ_' num2str(ID,'%02u') '.png'],'-dpng','-r600')
    
end

function [col] = color_error(di)

col = zeros(length(di),3);
for k = 1:length(di)
    col(k,:) = col_di(di(k));
end

end
function [col] = col_di(di)
% k = di/100;    
% col = [1 1-k 0];
if di==1
    col = 0.3*[1 1 1];
else
    col = 0.7*[1 1 1];
end
end
function [lev] = risklevel(epg)
lev = 0.*epg + 5;
lev(epg<200) = 4;
lev(epg<100) = 3;
lev(epg<70) = 2;
lev(epg<30) = 1;
end
function [ovl] = overlap(x,ix,y,iy)
ovl = (x+ix > y-iy) & (x-ix < y+iy);
end
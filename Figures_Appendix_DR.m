close all
clear all

close all
clear all

load('./Data/tables.mat')

eq = readtable('./Data/relacio_aga.xlsx');
prop = xlsread('./Figures&Tables/Table_incidence_T2.xlsx','Numbers','L:L');
prop = [prop(2:end); prop(1)];
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

A14_PCR = zeros(length(time),length(u)+1);
A14_EMR = zeros(length(time),length(u)+1);
iA14_PCR = zeros(length(time),length(u)+1);
iA14_EMR = zeros(length(time),length(u)+1);

cas_PCR = zeros(length(time),length(u)+1);
cas_EMR = zeros(length(time),length(u)+1);
icas_PCR = zeros(length(time),length(u)+1);
icas_EMR = zeros(length(time),length(u)+1);

cas7_PCR = zeros(length(time),length(u)+1);
cas7_EMR = zeros(length(time),length(u)+1);
icas7_PCR = zeros(length(time),length(u)+1);
icas7_EMR = zeros(length(time),length(u)+1);

rho_PCR = zeros(length(time),length(u)+1);
rho_EMR = zeros(length(time),length(u)+1);
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
        
        [phat, pci] = binofit(cas_PCR(k,i),pop(i),0.05);
        icas_PCR(k,i) = mean(abs(pci-phat)*pop(i));
        [phat, pci] = binofit(cas_EMR(k,i),pop(i),0.05);
        icas_EMR(k,i) = mean(abs(pci-phat)*pop(i));
        
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

u = [u;0];
Name = cell(length(u),1);
for k = 1:length(u)-1
    Name(k) = eq.aga(eq.id_aga==u(k));
end
Name{end} = 'CATALUNYA';

tt_Jul = datenum('01/07/2020','dd/mm/yyyy'):datenum('31/07/2020','dd/mm/yyyy');
id_t = ismember(time,tt_Jul);
A14_PCR_Jul = A14_PCR(id_t,:);
A14_EMR_Jul = A14_EMR(id_t,:);
rho_PCR_Jul = rho_PCR(id_t,:);
rho_EMR_Jul = rho_EMR(id_t,:);
iA14_PCR_Jul = iA14_PCR(id_t,:);
iA14_EMR_Jul = iA14_EMR(id_t,:);
irho_PCR_Jul = irho_PCR(id_t,:);
irho_EMR_Jul = irho_EMR(id_t,:);

tt_Aug = datenum('01/08/2020','dd/mm/yyyy'):datenum('31/08/2020','dd/mm/yyyy');
id_t = ismember(time,tt_Aug);
A14_PCR_Aug = A14_PCR(id_t,:);
A14_EMR_Aug = A14_EMR(id_t,:);
rho_PCR_Aug = rho_PCR(id_t,:);
rho_EMR_Aug = rho_EMR(id_t,:);
iA14_PCR_Aug = iA14_PCR(id_t,:);
iA14_EMR_Aug = iA14_EMR(id_t,:);
irho_PCR_Aug = irho_PCR(id_t,:);
irho_EMR_Aug = irho_EMR(id_t,:);

AGA = {};
id_aga = [u;0];
eq = readtable('./Data/relacio_aga.xlsx');
for k = 1:length(u)
    AGA = [AGA; eq.aga(eq.id_aga==u(k))];
end
AGA = [AGA; {'CATALUNYA'}];

%%
LW = 1;

for id = 1:44
    
    f = figure(1);
    
    f.Position = [33   180   880   360];
    col = [0 0.447 0.741; 0.5 0.17 0.05];
    
    clf
    
    rx = 0.03;
    lx = 0.08;
    mx = 0.08;
    uy = 0.15;
    dy = 0.2;
    DX = (1-rx-mx-lx)/2;
    DY = 1-uy-dy;
    
    errorbar(A14_PCR_Jul(:,id),rho_PCR_Jul(:,id),irho_PCR_Jul(:,id),irho_PCR_Jul(:,id),iA14_PCR_Jul(:,id),iA14_PCR_Jul(:,id),'o','Color',col(1,:),'linewidth',LW)
    hold on
    errorbar(A14_EMR_Jul(:,id)*prop(id),rho_EMR_Jul(:,id),irho_EMR_Jul(:,id),irho_EMR_Jul(:,id),iA14_EMR_Jul(:,id)*prop(id),iA14_EMR_Jul(:,id)*prop(id),'o','Color',col(2,:),'linewidth',LW)
    ax = gca;
    supxJ = ax.XLim(2);
    
    clf
    
    errorbar(A14_PCR_Aug(:,id),rho_PCR_Aug(:,id),irho_PCR_Aug(:,id),irho_PCR_Aug(:,id),iA14_PCR_Aug(:,id),iA14_PCR_Aug(:,id),'o','Color',col(1,:),'linewidth',LW)
    hold on
    errorbar(A14_EMR_Aug(:,id)*prop(id),rho_EMR_Aug(:,id),irho_EMR_Aug(:,id),irho_EMR_Aug(:,id),iA14_EMR_Aug(:,id)*prop(id),iA14_EMR_Aug(:,id)*prop(id),'o','Color',col(2,:),'linewidth',LW)
    ax = gca;
    supxA = ax.XLim(2);
    
    clf
    
    ax = axes('Units','normalized','Position',[lx dy DX DY]);
    
    plot(A14_PCR_Jul(:,id),rho_PCR_Jul(:,id),'-o','Color',col(1,:),'linewidth',LW)
    hold on
    plot(prop(id)*A14_EMR_Jul(:,id),rho_EMR_Jul(:,id),'-o','Color',col(2,:),'linewidth',LW)
    fons(supxJ)
    
    n = size(A14_EMR_Jul,1);
    for i = 1:n
        quadrats(A14_PCR_Jul(i,id),iA14_PCR_Jul(i,id),rho_PCR_Jul(i,id),irho_PCR_Jul(i,id),col(1,:),0.5,0.5)
    end
    for i = 1:n
        quadrats(A14_EMR_Jul(i,id)*prop(id),iA14_EMR_Jul(i,id)*prop(id),rho_EMR_Jul(i,id),irho_EMR_Jul(i,id),col(2,:),0.5,0.5)
    end
    
    plot([0 supxJ],[1 1],'--k')
    
    plot(A14_PCR_Jul(:,id),rho_PCR_Jul(:,id),'-o','Color',col(1,:),'linewidth',LW)
    plot(prop(id)*A14_EMR_Jul(:,id),rho_EMR_Jul(:,id),'-o','Color',col(2,:),'linewidth',LW)
    
    errorbar(A14_PCR_Jul(end,id),rho_PCR_Jul(end,id),irho_PCR_Jul(end,id),irho_PCR_Jul(end,id),iA14_PCR_Jul(end,id),iA14_PCR_Jul(end,id),'o','Color',col(1,:),'linewidth',LW)
    errorbar(A14_EMR_Jul(end,id)*prop(id),rho_EMR_Jul(end,id),irho_EMR_Jul(end,id),irho_EMR_Jul(end,id),iA14_EMR_Jul(end,id)*prop(id),iA14_EMR_Jul(end,id)*prop(id),'o','Color',col(2,:),'linewidth',LW)
    
    plot(A14_PCR_Jul(end,id),rho_PCR_Jul(end,id),'-o','Color',col(1,:),'MarkerFaceColor',0.2*col(1,:)+0.8,'linewidth',LW)
    plot(A14_PCR_Jul(1,id),rho_PCR_Jul(1,id),'-o','Color',col(1,:),'MarkerFaceColor',0.5*col(1,:),'linewidth',LW)
    plot(prop(id)*A14_EMR_Jul(end,id),rho_EMR_Jul(end,id),'-o','Color',col(2,:),'MarkerFaceColor',0.2*col(2,:)+0.8,'linewidth',LW)
    plot(prop(id)*A14_EMR_Jul(1,id),rho_EMR_Jul(1,id),'-o','Color',col(2,:),'MarkerFaceColor',0.5*col(2,:),'linewidth',LW)
    
    ax0 = gca;
    ax0.LineWidth = 1.5;
    ax0.Box = 'off';
    ax0.TickDir = 'out';
    axis([0 supxJ 0 4])
    xlabel('Active cases per 10^5 inh. (A_1_4)')
    ylabel('Empiric propagation')
    
    xl = xlim;
    yl = ylim;
    plot(xl,yl(1)*[1 1],'-k','linewidth',1)
    plot(xl(1)*[1 1],yl,'-k','linewidth',1)
    
    l1A=legend('PCR','EMR','Location','northeast');
    %l1A.Position(1:2) = [ax.Position(1)+0.02, ax.Position(2)+ax.Position(4)-0.02-l1A.Position(4)];
    text(0.03,1,'A','Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',24)
    
    ax1 = axes('Units','normalized','Position',[-1 -1 0.5 0.5]);
    plot(1,1,'ok','MarkerFaceColor',0.25*[1 1 1])
    hold on
    plot(1,1,'ok','MarkerFaceColor',0.9*[1 1 1])
    l2A = legend(datestr(tt_Jul(1)),datestr(tt_Jul(end)));
    l2A.Position(1:3) = [ax.Position(1), ax.Position(2)+ax.Position(4)+0.02 ,0.1336];
    l2A.LineWidth = 1.5;
    
    ax = axes('Units','normalized','Position',[lx+mx+DX dy DX DY]);
    
    plot(A14_PCR_Aug(:,id),rho_PCR_Aug(:,id),'-o','Color',col(1,:),'linewidth',LW)
    hold on
    plot(prop(id)*A14_EMR_Aug(:,id),rho_EMR_Aug(:,id),'-o','Color',col(2,:),'linewidth',LW)
    fons(supxA)
    
    n = size(A14_EMR_Aug,1);
    for i = 1:n
        quadrats(A14_PCR_Aug(i,id),iA14_PCR_Aug(i,id),rho_PCR_Aug(i,id),irho_PCR_Aug(i,id),col(1,:),0.5,0.5)
    end
    for i = 1:n
        quadrats(A14_EMR_Aug(i,id)*prop(id),iA14_EMR_Aug(i,id)*prop(id),rho_EMR_Aug(i,id),irho_EMR_Aug(i,id),col(2,:),0.5,0.5)
    end
    
    plot([0 supxA],[1 1],'--k')
    
    plot(A14_PCR_Aug(:,id),rho_PCR_Aug(:,id),'-o','Color',col(1,:),'linewidth',LW)
    plot(prop(id)*A14_EMR_Aug(:,id),rho_EMR_Aug(:,id),'-o','Color',col(2,:),'linewidth',LW)
    
    errorbar(A14_PCR_Aug(end,id),rho_PCR_Aug(end,id),irho_PCR_Aug(end,id),irho_PCR_Aug(end,id),iA14_PCR_Aug(end,id),iA14_PCR_Aug(end,id),'o','Color',col(1,:),'linewidth',LW)
    errorbar(A14_EMR_Aug(end,id)*prop(id),rho_EMR_Aug(end,id),irho_EMR_Aug(end,id),irho_EMR_Aug(end,id),iA14_EMR_Aug(end,id)*prop(id),iA14_EMR_Aug(end,id)*prop(id),'o','Color',col(2,:),'linewidth',LW)
    
    plot(A14_PCR_Aug(end,id),rho_PCR_Aug(end,id),'-o','Color',col(1,:),'MarkerFaceColor',0.2*col(1,:)+0.8,'linewidth',LW)
    plot(A14_PCR_Aug(1,id),rho_PCR_Aug(1,id),'-o','Color',col(1,:),'MarkerFaceColor',0.5*col(1,:),'linewidth',LW)
    plot(prop(id)*A14_EMR_Aug(end,id),rho_EMR_Aug(end,id),'-o','Color',col(2,:),'MarkerFaceColor',0.2*col(2,:)+0.8,'linewidth',LW)
    plot(prop(id)*A14_EMR_Aug(1,id),rho_EMR_Aug(1,id),'-o','Color',col(2,:),'MarkerFaceColor',0.5*col(2,:),'linewidth',LW)
    
    ax0 = gca;
    ax0.LineWidth = 1.5;
    ax0.Box = 'off';
    ax0.TickDir = 'out';
    axis([0 supxA 0 4])
    xlabel('Active cases per 10^5 inh. (A_1_4)')
    ylabel('Empiric propagation')
    
    xl = xlim;
    yl = ylim;
    plot(xl,yl(1)*[1 1],'-k','linewidth',1)
    plot(xl(1)*[1 1],yl,'-k','linewidth',1)
    
    l1A=legend('PCR','EMR','Location','northeast');
    %l1A.Position(1:2) = [ax.Position(1)+ax.Position(3)-l1A.Position(3)-0.02, ax.Position(2)+ax.Position(4)-0.02-l1A.Position(4)];
    text(0.03,1,'B','Units','normalized','HorizontalAlignment','left','VerticalAlignment','top','FontSize',24)
    
    ax1 = axes('Units','normalized','Position',[-1 -1 0.5 0.5]);
    plot(1,1,'ok','MarkerFaceColor',0.25*[1 1 1])
    hold on
    plot(1,1,'ok','MarkerFaceColor',0.9*[1 1 1])
    l2A = legend(datestr(tt_Aug(1)),datestr(tt_Aug(end)));
    l2A.Position(1:3) = [ax.Position(1)+ax.Position(3)-0.1336, ax.Position(2)+ax.Position(4)+0.02 ,0.1336];
    l2A.LineWidth = 1.5;
    
    a = annotation('textbox',[lx+DX 0.87 0.1 0.1],'String',[AGA{id} ', population: ' num2str(pop(id))],'FitBoxToText','on','EdgeColor','none','FontSize',15);
    drawnow;
    a.Position(1) = 0.53- a.Position(3)/2;
    drawnow;
    
    print(f,['./Figures&Tables/Appendix_DR/DR_' num2str(id,'%02u') '.png'],'-dpng','-r600')
    
end

function [] = fons(supx)
rh = 0:0.001:3;
ar = 0:0.01:(supx+2);
[RH,AR] = meshgrid(rh,ar);
EPG = RH.*AR;

EPG(EPG>30 & EPG<40) =40;
EPG(EPG>100)=100;
x = 0:0.001:1;
x=(sqrt(x))';
a=0:(supx+20);
b=(0:5);
cmap = flip([0.*x+1 x 0.*x; flip(x) 0.*x+1 0.*x]);
ax=gca;
ax.YDir='normal';
im=imagesc(a,b,EPG);
ax.YDir='normal';
colormap(cmap)
im.AlphaData = 0.4;
end
function [] = quadrats(x,dx,y,dy,col,wt,al)

N = length(x);

for k = 1:N
    
    fill([x(k)-dx(k) x(k)-dx(k) x(k)+dx(k) x(k)+dx(k) x(k)-dx(k)],...
        [y(k)-dy(k) y(k)+dy(k) y(k)+dy(k) y(k)-dy(k) y(k)-dy(k)],...
        (1-wt)*col+wt,'EdgeColor','none','FaceAlpha',al)
    
end

end

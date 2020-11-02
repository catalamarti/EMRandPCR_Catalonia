close all
clear all

load('./Data/tables.mat')

tt = datenum('27/08/2020','dd/mm/yyyy'):datenum('13/09/2020','dd/mm/yyyy');


u = unique(T0827.id_AGA);
pop = 0.*u;
for i = 1:length(u)
    pop(i) = unique(T0827.popAGA(T0827.id_AGA==u(i)));
end
pop = [sum(pop); pop];

for k = 1:length(tt)
    eval(['T = T' datestr(tt(k),'mmdd') ';']);
    dies = unique(datenum(T.data));
    PCR = 0.*dies;
    EMR = 0.*dies;
    for i = 1:length(dies)
        PCR(i) = sum(T.PCRcas(datenum(T.data)==dies(i)));
        EMR(i) = sum(T.ClinicsCovid(datenum(T.data)==dies(i)));
    end
    eval(['t' datestr(tt(k),'mmdd') ' = table(dies,PCR,EMR);']);
end

delay_PCR = [];
delay_EMR = [];

for k = 2:length(tt)
    
    eval(['Told = t' datestr(tt(k)-1,'mmdd') ';']);
    eval(['Tnow = t' datestr(tt(k),'mmdd') ';']);
    
    nous = abs(Tnow.PCR - [Told.PCR;0]);
    dies_retard = tt(k)-Tnow.dies;
    N = sum(nous);
    ret = zeros(N,1);
    j = 0;
    for i = 1:length(nous)
        ret((1:nous(i))+j) = dies_retard(i);
        j = j + nous(i);
    end
    delay_PCR = [delay_PCR; ret];
    
    nous = abs(Tnow.EMR - [Told.EMR;0]);
    dies_retard = tt(k)-Tnow.dies;
    N = sum(nous);
    ret = zeros(N,1);
    j = 0;
    for i = 1:length(nous)
        ret((1:nous(i))+j) = dies_retard(i);
        j = j + nous(i);
    end
    delay_EMR = [delay_EMR; ret];
    
end

delay_PCR = delay_PCR(delay_PCR<30);
delay_EMR = delay_EMR(delay_EMR<30);

ncPCR = hist(delay_PCR(delay_PCR<30),0:200);
ncPCR = ncPCR/sum(ncPCR);
ncEMR = hist(delay_EMR(delay_EMR<30),0:200);
ncEMR = ncEMR/sum(ncEMR);

col = lines(4);

f = figure(1);
clf

plot(0:200,ncPCR,'Color',col(1,:),'linewidth',2)
hold on
plot(0:200,ncEMR,'Color',col(2,:),'linewidth',2)

plot(mean(delay_PCR(delay_PCR<30))*[1 1],[0 1],'Color',col(1,:),'linewidth',1)
plot(mean(delay_EMR(delay_EMR<30))*[1 1],[0 1],'Color',col(2,:),'linewidth',1)

ylabel('Fraction of reported cases')
xlabel('Reporting back day')

axis([0 10 0 1])

l=legend('PCR','EMR');
l.FontSize=12;

ax = gca;
ax.Box = 'off';
ax.TickDir = 'out';
ax.LineWidth = 1.5;

print(f,'./Figures&Tables/Figure_distributions_SF1.png','-dpng','-r600')

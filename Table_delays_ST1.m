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
    
    eval(['Tantic = t' datestr(tt(k)-1,'mmdd') ';']);
    eval(['Tnow = t' datestr(tt(k),'mmdd') ';']);
    
    nous = abs(Tnow.PCR - [Tantic.PCR;0]);
    dies_retard = tt(k)-Tnow.dies;
    N = sum(nous);
    ret = zeros(N,1);
    j = 0;
    for i = 1:length(nous)
        ret((1:nous(i))+j) = dies_retard(i);
        j = j + nous(i);
    end
    delay_PCR = [delay_PCR; ret];
    
    nous = abs(Tnow.EMR - [Tantic.EMR;0]);
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

Deleyas = [mean(delay_PCR) median(delay_PCR) find(cumsum(ncPCR)>0.95,1,'first') ...
    mean(delay_EMR) median(delay_EMR) find(cumsum(ncEMR)>0.95,1,'first')];

NcasosPCR = [length(delay_PCR); 0.*u];
NcasosEMR = [length(delay_EMR); 0.*u];

for ii = 1:length(u)
    
    for k = 1:length(tt)
        eval(['T = T' datestr(tt(k),'mmdd') ';']);
        dies = unique(datenum(T.data));
        PCR = 0.*dies;
        EMR = 0.*dies;
        for i = 1:length(dies)
            PCR(i) = sum(T.PCRcas(datenum(T.data)==dies(i) & T.id_AGA==u(ii)));
            EMR(i) = sum(T.ClinicsCovid(datenum(T.data)==dies(i) & T.id_AGA==u(ii)));
        end
        eval(['t' datestr(tt(k),'mmdd') ' = table(dies,PCR,EMR);']);
    end
    
    delay_PCR = [];
    delay_EMR = [];
    
    for k = 2:length(tt)
        
        eval(['Tantic = t' datestr(tt(k)-1,'mmdd') ';']);
        eval(['Tnow = t' datestr(tt(k),'mmdd') ';']);
        
        nous = abs(Tnow.PCR - [Tantic.PCR;0]);
        dies_retard = tt(k)-Tnow.dies;
        N = sum(nous);
        ret = zeros(N,1);
        j = 0;
        for i = 1:length(nous)
            ret((1:nous(i))+j) = dies_retard(i);
            j = j + nous(i);
        end
        delay_PCR = [delay_PCR; ret];
        
        nous = abs(Tnow.EMR - [Tantic.EMR;0]);
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
    
    NcasosPCR(ii+1) = length(delay_PCR);
    NcasosEMR(ii+1) = length(delay_EMR);
    
    ncPCR = hist(delay_PCR(delay_PCR<30),0:200);
    ncPCR = ncPCR/sum(ncPCR);
    ncEMR = hist(delay_EMR(delay_EMR<30),0:200);
    ncEMR = ncEMR/sum(ncEMR);
    
    Deleyas = [Deleyas; mean(delay_PCR) median(delay_PCR) find(cumsum(ncPCR)>0.95,1,'first') ...
    mean(delay_EMR) median(delay_EMR) find(cumsum(ncEMR)>0.95,1,'first')];
    
end

AGA = {'CATALUNYA'};
id_aga = [0;u];

eq = readtable('./Data/relacio_aga.xlsx');
for k = 1:length(u)
    AGA = [AGA; eq.aga(eq.id_aga==u(k))];
end

taula = table(AGA,Deleyas(:,1),Deleyas(:,2),Deleyas(:,3),Deleyas(:,4),Deleyas(:,5),Deleyas(:,6));
taula.Properties.VariableNames = {'Health_region','mean_PCR','median_PCR','95%_PCR','mean_EMR','median_EMR','95%_EMR'};
writetable(taula,'./Figures&Tables/Table_delays_SF1.xlsx')
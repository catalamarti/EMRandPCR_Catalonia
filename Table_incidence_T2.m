close all
clear all

load('./Data/tables.mat')


time = datenum('01/07/2020','dd/mm/yyyy'):datenum('31/08/2020','dd/mm/yyyy');
PCR = 0.*time;
EMR = 0.*time;
T = T0913;

u = unique(T.id_AGA(~isnan(T.id_AGA)));
pop = 0.*u;
for k = 1:length(u)
    pop(k) = unique(T.popAGA(T.id_AGA==u(k)));
end
pop = [sum(pop); pop];

inf_PCR = 0.*time;
sup_PCR = 0.*time;
inf_EMR = 0.*time;
sup_EMR = 0.*time;

for k = 1:length(time)
    PCR(k)      = sum(T.PCRcas(datenum(T.data)==time(k)));
    [ph,pci] = binofit(PCR(k),pop(1),0.05);
    inf_PCR(k) = (ph-pci(1))*pop(1);
    sup_PCR(k) = (pci(2)-ph)*pop(1);
    EMR(k) = sum(T.ClinicsCovid(datenum(T.data)==time(k)));
    inf_EMR(k) = (ph-pci(1))*pop(1);
    sup_EMR(k) = (pci(2)-ph)*pop(1);
end
s_inf_PCR = sqrt(sum(inf_PCR.^2));
s_inf_EMR = sqrt(sum(inf_EMR.^2));
s_sup_PCR = sqrt(sum(sup_PCR.^2));
s_sup_EMR = sqrt(sum(sup_EMR.^2));
p = sum(EMR)/sum(PCR);
sig_inf_p = 1/sum(PCR)*sqrt(p*p*s_inf_PCR^2+s_inf_EMR^2);
sig_sup_p = 1/sum(PCR)*sqrt(p*p*s_inf_PCR^2+s_sup_EMR^2);

Data = [sum(PCR) sum(PCR)/pop(1)*1e5 (sum(PCR)-s_inf_PCR)/pop(1)*1e5 ...
    (sum(PCR)+s_sup_PCR)/pop(1)*1e5 sum(EMR) sum(EMR)/pop(1)*1e5 (sum(EMR)-s_inf_EMR)/pop(1)*1e5 ...
    (sum(EMR)+s_sup_EMR)/pop(1)*1e5 p p-sig_inf_p p+sig_inf_p];

for i = 1:length(u)
    
    for k = 1:length(time)
        PCR(k)      = sum(T.PCRcas(T.id_AGA==u(i) & datenum(T.data)==time(k)));
        [ph,pci] = binofit(PCR(k),pop(i+1),0.05);
        inf_PCR(k) = (ph-pci(1))*pop(i+1);
        sup_PCR(k) = (pci(2)-ph)*pop(i+1);
        EMR(k) = sum(T.ClinicsCovid(T.id_AGA==u(i) & datenum(T.data)==time(k)));
        inf_EMR(k) = (ph-pci(1))*pop(i+1);
        sup_EMR(k) = (pci(2)-ph)*pop(i+1);
    end
    
    s_inf_PCR = sqrt(sum(inf_PCR.^2));
    s_inf_EMR = sqrt(sum(inf_EMR.^2));
    s_sup_PC = sqrt(sum(sup_PCR.^2));
    s_sup_EMR = sqrt(sum(sup_EMR.^2));
    p = sum(EMR)/sum(PCR);
    sig_inf_p = 1/sum(PCR)*sqrt(p*p*s_inf_PCR^2+s_inf_EMR^2);
    sig_sup_p = 1/sum(PCR)*sqrt(p*p*s_inf_PCR^2+s_sup_EMR^2);

    Data = [Data;sum(PCR) sum(PCR)/pop(i+1)*1e5 (sum(PCR)-s_inf_PCR)/pop(i+1)*1e5 ...
    (sum(PCR)+s_sup_PCR)/pop(i+1)*1e5 sum(EMR) sum(EMR)/pop(i+1)*1e5 (sum(EMR)-s_inf_EMR)/pop(i+1)*1e5 ...
    (sum(EMR)+s_sup_EMR)/pop(i+1)*1e5 p p-sig_inf_p p+sig_inf_p];
    
end

AGA = {'CATALUNYA'};
id_aga = [0;u];

eq = readtable('./Data/relacio_aga.xlsx');
for k = 1:length(u)
    AGA = [AGA; eq.aga(eq.id_aga==u(k))];
end

taula = table(AGA,id_aga,pop,Data(:,1),Data(:,2),Data(:,3),Data(:,4),Data(:,5),Data(:,6),Data(:,7),Data(:,8),Data(:,9),Data(:,10),Data(:,11));
taula.Properties.VariableNames = {'Health_region','id','Population','PCR','PCR_CI',...
    'PCR_CI_inf','PCR_CI_sup','EMR','EMR_CI','EMR_CI_inf','EMR_CI_sup','EMR/PCR',...
    'EMR/PCR_inf','EMR/PCR_sup'};
writetable(taula,'./Figures&Tables/Table_incidence_T2.xlsx','sheet','Numbers')

PCR_int = cell(length(AGA),1);
EMR_int = cell(length(AGA),1);
quotient_int = cell(length(AGA),1);

for k = 1:length(AGA)
    form = '%5.1f';
    PCR_int{k} = [num2str(Data(k,3),form) ' - ' num2str(Data(k,4),form)];
    EMR_int{k} = [num2str(Data(k,7),form) ' - ' num2str(Data(k,8),form)];
    form = '%4.2f';
    quotient_int{k} = [num2str(Data(k,10),form) ' - ' num2str(Data(k,11),form)];
end

taula = table(AGA,pop,Data(:,1),Data(:,2),PCR_int,Data(:,5),Data(:,6),EMR_int,Data(:,9),quotient_int);
taula.Properties.VariableNames = {'Health_region','Population','PCR','PCR_CI',...
    'PCR_CI_int','EMR','EMR_CI','EMR_CI_int','EMR/PCR','EMR/PCR_int'};
writetable(taula,'./Figures&Tables/Table_incidence_T2.xlsx','sheet','Table')

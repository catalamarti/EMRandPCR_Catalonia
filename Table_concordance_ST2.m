close all
clear all

load('./Data/tables.mat')
eq = readtable('./Data/relacio_aga.xlsx');
prop = xlsread('./Figures&Tables/Table_incidence_T2.xlsx','Numbers','L:L');
prop = [prop(2:end);prop(1)];
prop = 1./prop;

time = datenum('01/06/2020','dd/mm/yyyy'):datenum('31/08/2020','dd/mm/yyyy');
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

EPG_CC = (prop.*EPG_CC')';
iEPG_CC = (prop.*iEPG_CC')';

time(1:30)=[];
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
            eval([s1 s2 '_' s3 '(1:30,:) = [];'])
        end
    end
end

AGA = {'CATALUNYA'};
id_aga = [0;u];
for k = 1:length(u)
    AGA = [AGA; eq.aga(eq.id_aga==u(k))];
end

A = zeros(44,18);
for j = 0:43
    if j == 0
        idd = 44;
    else
        idd = j;
    end
    epg_0 = EPG_Pc(:,idd);
    epg_i = EPG_CC(:,idd);
    epg_0(epg_0>200) = -5;
    epg_0(epg_0>100) = -4;
    epg_0(epg_0>70) = -3;
    epg_0(epg_0>30) = -2;
    epg_0(epg_0>=0) = -1;
    epg_0 = abs(epg_0);
    epg_i(epg_i>200) = -5;
    epg_i(epg_i>100) = -4;
    epg_i(epg_i>70) = -3;
    epg_i(epg_i>30) = -2;
    epg_i(epg_i>=0) = -1;
    epg_i = abs(epg_i);
    
    for i = 0:5
        
        if i == 0
            id = true(62,1);
        else
            eval(['id = epg_0==' num2str(i) ';']);
        end
        n = 0;
        EPG_0 = epg_0(id);
        EPG_i = epg_i(id);
        N = length(EPG_0);
        epg_val_0 = EPG_Pc(id,idd);
        epg_val_i = EPG_CC(id,idd);
        iepg_val_0 = iEPG_Pc(id,idd);
        iepg_val_i = iEPG_CC(id,idd);
        for k = 1:N
            if EPG_0(k)==EPG_i(k)
                n = n+1;
            elseif abs(epg_val_0(k)-epg_val_i(k))<=iepg_val_0(k)+iepg_val_i(k)
                n = n+1;
                %fprintf('PCR:%4.1f+-%3.1f vs CC:%4.1f+-%3.1f\n',epg_val_0(k),iepg_val_0(k),epg_val_i(k),iepg_val_i(k))
            end
        end
        A(j+1,3*(i+1)-2) = N;
        A(j+1,3*(i+1)-1) = n;
        A(j+1,3*(i+1)) = n/N;
        
    end
    
end

xlswrite('./Figures&Tables/Table_concordance_ST2.xlsx',AGA,'Numbers','A2')
xlswrite('./Figures&Tables/Table_concordance_ST2.xlsx',pop,'Numbers','B3')
xlswrite('./Figures&Tables/Table_concordance_ST2.xlsx',A,'Numbers','C2')

close all
clear all

f = dir('./Figures&Tables/Appendix_DR/DR_*');

for k = 1:round(length(f)/4)
    
    I1 = imread([f(4*(k-1)+1).folder '/' f(4*(k-1)+1).name]);
    I2 = imread([f(4*(k-1)+2).folder '/' f(4*(k-1)+2).name]);
    I3 = imread([f(4*(k-1)+3).folder '/' f(4*(k-1)+3).name]);
    I4 = imread([f(4*(k-1)+4).folder '/' f(4*(k-1)+4).name]);
    
    h = 2;
    lx = 3/h;
    uy = 2.7/h;
    DX = 15/h;
    DY = 6/h;
    
    fi = figure(1);
    clf
    fi.Units = 'centimeters';
    fi.Position = [0 0 21/h 29.7/h];
    ax = axes('Units','centimeters','Position',[lx uy+3*DY DX DY]);
    imagesc(I1)
    axis off
    ax = axes('Units','centimeters','Position',[lx uy+2*DY DX DY]);
    imagesc(I2)
    axis off
    ax = axes('Units','centimeters','Position',[lx uy+1*DY DX DY]);
    imagesc(I3)
    axis off
    ax = axes('Units','centimeters','Position',[lx uy+0*DY DX DY]);
    imagesc(I4)
    axis off
    
    print(fi,['./Figures&Tables/Appendix_DR/Fig' num2str(k,'%02u') '.pdf'],'-dpdf','-fillpage')
    
end

append_pdfs('./Figures&Tables/Appendix_DR.pdf','./Figures&Tables/Appendix_DR/Fig01.pdf',...
    './Figures&Tables/Appendix_DR/Fig02.pdf','./Figures&Tables/Appendix_DR/Fig03.pdf',...
    './Figures&Tables/Appendix_DR/Fig04.pdf','./Figures&Tables/Appendix_DR/Fig05.pdf',...
    './Figures&Tables/Appendix_DR/Fig06.pdf','./Figures&Tables/Appendix_DR/Fig07.pdf',...
    './Figures&Tables/Appendix_DR/Fig08.pdf','./Figures&Tables/Appendix_DR/Fig09.pdf',...
    './Figures&Tables/Appendix_DR/Fig10.pdf','./Figures&Tables/Appendix_DR/Fig11.pdf')

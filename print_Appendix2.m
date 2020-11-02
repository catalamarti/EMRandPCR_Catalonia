close all
clear all

load('./Data/names.mat')
f = dir('./Figures&Tables/Appendix2/*.png');

for k = 1:22
    
    i1 = 2*(k-1)+1;
    i2 = 2*k;
    F1 = imread([f(i1).folder '/' f(i1).name]);
    if k < 22
        F2 = imread([f(i2).folder '/' f(i2).name]);
    end
    
    h=figure(1);
    clf
    K = 2;
    set(h,'Units','centimeters');
    set(h,'Position',[0 0 21/K 29.7/K])
    set(h,'PaperPositionMode','Auto')
    set(h,'PaperSize',[21 29.7])
    
    dy = 4.2;
    my = 2.2;
    lx = 3;
    
    LY = 10.3;
    LX = 15;
    TH = 25;
    
    axes('Units','centimeters','Position',[lx dy+my+LY LX LY]/K)
    imagesc(F1)
    axis equal
    axis off
    text(0,-0.05,['{\bfFigure ' num2str(i1) '}. Weekly average of daily new cases (A), \rho_7 empirical propagation (B), A_1_4 attack'],'units','normalized')
    if length(AGA{i1})>TH
        id = strfind(AGA{i1},' ');
        ii = id(find(id<TH,1,'last'));
        nam1 = AGA{i1}(1:ii);
        nam2 = AGA{i1}(ii+1:end);
        text(0,-0.10,['rate (C) and EPG empirical propagation growth (D) over time for ' nam1],'units','normalized')
        text(0,-0.15,[nam2 '.'],'units','normalized')
    else
        text(0,-0.10,['rate (C) and EPG empirical propagation growth (D) over time for ' AGA{i1} '.'],'units','normalized')
    end
    if k<22
        axes('Units','centimeters','Position',[lx dy LX LY]/K)
        imagesc(F2)
        axis equal
        axis off
        text(0,-0.05,['{\bfFigure ' num2str(i2) '}. Weekly average of daily new cases (A), \rho_7 empirical propagation (B), A_1_4 attack'],'units','normalized')
        if length(AGA{i2})>TH
            id = strfind(AGA{i2},' ');
            ii = id(find(id<TH,1,'last'));
            nam1 = AGA{i2}(1:ii);
            nam2 = AGA{i2}(ii+1:end);
            text(0,-0.10,['rate (C) and EPG empirical propagation growth (D) over time for ' nam1],'units','normalized')
            text(0,-0.15,[nam2 '.'],'units','normalized')
        else
            text(0,-0.10,['rate (C) and EPG empirical propagation growth (D) over time for ' AGA{i2} '.'],'units','normalized')
            
        end
    end
    
    print(h,['./Figures&Tables/Appendix2/Fig' num2str(k,'%02u') '.pdf'],'-fillpage','-dpdf','-r0')
    
end

append_pdfs('./Figures&Tables/Appendix2.pdf','./Figures&Tables/Appendix2/Fig01.pdf',...
    './Figures&Tables/Appendix2/Fig02.pdf','./Figures&Tables/Appendix2/Fig03.pdf',...
    './Figures&Tables/Appendix2/Fig04.pdf','./Figures&Tables/Appendix2/Fig05.pdf',...
    './Figures&Tables/Appendix2/Fig06.pdf','./Figures&Tables/Appendix2/Fig07.pdf',...
    './Figures&Tables/Appendix2/Fig08.pdf','./Figures&Tables/Appendix2/Fig09.pdf',...
    './Figures&Tables/Appendix2/Fig10.pdf','./Figures&Tables/Appendix2/Fig11.pdf',...
    './Figures&Tables/Appendix2/Fig12.pdf','./Figures&Tables/Appendix2/Fig13.pdf',...
    './Figures&Tables/Appendix2/Fig14.pdf','./Figures&Tables/Appendix2/Fig15.pdf',...
    './Figures&Tables/Appendix2/Fig16.pdf','./Figures&Tables/Appendix2/Fig17.pdf',...
    './Figures&Tables/Appendix2/Fig18.pdf','./Figures&Tables/Appendix2/Fig19.pdf',...
    './Figures&Tables/Appendix2/Fig20.pdf','./Figures&Tables/Appendix2/Fig21.pdf',...
    './Figures&Tables/Appendix2/Fig22.pdf')


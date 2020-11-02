close all
clear all

load('./Data/names.mat')
f = dir('./Figures&Tables/Appendix_EPG_scatter/*.png');

for k = 1:8
    
    i1 = 6*(k-1)+1;
    i2 = 6*(k-1)+2;
    i3 = 6*(k-1)+3;
    i4 = 6*(k-1)+4;
    i5 = 6*(k-1)+5;
    i6 = 6*k;
    F1 = imread([f(i1).folder '/' f(i1).name]); F1 = F1(:,400:3000,:);
    F2 = imread([f(i2).folder '/' f(i2).name]); F2 = F2(:,400:3000,:);
    if k < 8
        F3 = imread([f(i3).folder '/' f(i3).name]); F3 = F3(:,400:3000,:);
        F4 = imread([f(i4).folder '/' f(i4).name]); F4 = F4(:,400:3000,:);
        F5 = imread([f(i5).folder '/' f(i5).name]); F5 = F5(:,400:3000,:);
        F6 = imread([f(i6).folder '/' f(i6).name]); F6 = F6(:,400:3000,:);
    end
    
    h=figure(1);
    clf
    K = 2;
    set(h,'Units','centimeters');
    set(h,'Position',[0 0 21/K 29.7/K])
    set(h,'PaperPositionMode','Auto')
    set(h,'PaperSize',[21 29.7])
    
    uy = 0;
    dy = 2;
    my = 1;
    lx = 2;
    rx = 2;
    mx = 2;
    fs = 8;
    b = 0;
    
    LY = (29.7-uy-2*my-dy)/3;
    LX = (21-lx-rx-mx)/2;
    TH = 25;
    
    axes('Units','centimeters','Position',[lx dy+2*my+2*LY LX LY]/K)
    imagesc(F1)
    axis equal
    axis off
    caption(i1,b,fs,AGA)
    
    axes('Units','centimeters','Position',[lx+LX+mx dy+2*my+2*LY LX LY]/K)
    imagesc(F2)
    axis equal
    axis off
    caption(i2,b,fs,AGA)
    
    if k<8
        
        axes('Units','centimeters','Position',[lx dy+my+LY LX LY]/K)
        imagesc(F3)
        axis equal
        axis off
        caption(i3,b,fs,AGA)
        
        axes('Units','centimeters','Position',[lx+LX+mx dy+my+LY LX LY]/K)
        imagesc(F4)
        axis equal
        axis off
        caption(i4,b,fs,AGA)
        
        axes('Units','centimeters','Position',[lx dy LX LY]/K)
        imagesc(F5)
        axis equal
        axis off
        caption(i5,b,fs,AGA)
        
        axes('Units','centimeters','Position',[lx+LX+mx dy LX LY]/K)
        imagesc(F6)
        axis equal
        axis off
        caption(i6,b,fs,AGA)
        
    end
    
    print(h,['./Figures&Tables/Appendix_EPG_scatter/Fig' num2str(k,'%02u') '.pdf'],'-fillpage','-dpdf','-r0')
    
end

delete('./Figures&Tables/Appendix3.pdf')
append_pdfs('./Figures&Tables/Appendix3.pdf','./Figures&Tables/Appendix_EPG_scatter/Fig01.pdf',...
    './Figures&Tables/Appendix_EPG_scatter/Fig02.pdf','./Figures&Tables/Appendix_EPG_scatter/Fig03.pdf',...
    './Figures&Tables/Appendix_EPG_scatter/Fig04.pdf','./Figures&Tables/Appendix_EPG_scatter/Fig05.pdf',...
    './Figures&Tables/Appendix_EPG_scatter/Fig06.pdf','./Figures&Tables/Appendix_EPG_scatter/Fig07.pdf',...
    './Figures&Tables/Appendix_EPG_scatter/Fig08.pdf')

function [] = caption(i,b,fs,AGA)
TH = 15;
text(0,b,['{\bfFigure ' num2str(i) '}. PCR- and EMR-EPG risk. Coincidences are marked'],'units','normalized','FontSize',fs)
text(0,b-0.05,['with dark gray and not coincidences in light gray. EPG measured'],'units','normalized','FontSize',fs)
if length(AGA{i})>TH
    id = strfind(AGA{i},' ');
    ii = id(find(id<TH,1,'last'));
    nam1 = AGA{i}(1:ii);
    nam2 = AGA{i}(ii+1:end);
    text(0,b-0.05*2,['between July 1st and August 31st. EPG for CMA: ' nam1],'units','normalized','FontSize',fs)
    text(0,b-0.05*3,[nam2 '.'],'units','normalized','FontSize',fs)
else
    text(0,b-0.05*2,['between July 1st and August 31st. EPG for CMA: ' AGA{i} '.'],'units','normalized','FontSize',fs)
end
end
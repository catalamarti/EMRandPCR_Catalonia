close all
clear all

load('./Data/names.mat')
f = dir('./Figures&Tables/Appendix_DR/*.png');

for k = 1:11
    
    i1 = 4*(k-1)+1;
    i2 = 4*(k-1)+2;
    i3 = 4*(k-1)+3;
    i4 = 4*(k-1)+4;
    F1 = imread([f(i1).folder '/' f(i1).name]);
    F2 = imread([f(i2).folder '/' f(i2).name]);
    F3 = imread([f(i3).folder '/' f(i3).name]);
    F4 = imread([f(i4).folder '/' f(i4).name]);
    
    h=figure(1);
    clf
    K = 2;
    set(h,'Units','centimeters');
    set(h,'Position',[0 0 21/K 29.7/K])
    set(h,'PaperPositionMode','Auto')
    set(h,'PaperSize',[21 29.7])
    
    uy = 0;
    dy = 1;
    my = 1;
    lx = 2;
    rx = 2;
    mx = 2;
    fs = 10;
    b = 0.03;
    
    LY = (29.7-uy-3*my-dy)/4;
    LX = (21-lx-rx-mx);
    TH = 25;
    
    axes('Units','centimeters','Position',[lx dy+3*my+3*LY LX LY]/K)
    imagesc(F1)
    axis equal
    axis off
    caption(i1,b,fs,AGA)
    
    axes('Units','centimeters','Position',[lx dy+2*(my+LY) LX LY]/K)
    imagesc(F2)
    axis equal
    axis off
    caption(i2,b,fs,AGA)
    
    axes('Units','centimeters','Position',[lx dy+my+LY LX LY]/K)
    imagesc(F3)
    axis equal
    axis off
    caption(i3,b,fs,AGA)
    
    axes('Units','centimeters','Position',[lx dy LX LY]/K)
    imagesc(F4)
    axis equal
    axis off
    caption(i4,b,fs,AGA)
    
    print(h,['./Figures&Tables/Appendix_DR/Fig' num2str(k,'%02u') '.pdf'],'-fillpage','-dpdf','-r0')
    
end

delete('./Figures&Tables/Appendix4.pdf')
append_pdfs('./Figures&Tables/Appendix4.pdf','./Figures&Tables/Appendix_DR/Fig01.pdf',...
    './Figures&Tables/Appendix_DR/Fig02.pdf','./Figures&Tables/Appendix_DR/Fig03.pdf',...
    './Figures&Tables/Appendix_DR/Fig04.pdf','./Figures&Tables/Appendix_DR/Fig05.pdf',...
    './Figures&Tables/Appendix_DR/Fig06.pdf','./Figures&Tables/Appendix_DR/Fig07.pdf',...
    './Figures&Tables/Appendix_DR/Fig08.pdf')

function [] = caption(i,b,fs,AGA)
TH = 95;
str = ['{\bfFigure ' num2str(i) '}. Risk diagram for the evolution of the COVID-19 pandemic in ' AGA{i} ' based on EMR (red) and on PCR (blue) cases for the month of July (left) and August 2020 (right).'];
ii = strfind(str,' ');
ii = ii(find(ii<TH,1,'last'));
str1 = str(1:ii);
str2 = str(ii+1:end);
text(0,b,str1,'units','normalized','FontSize',fs)
if length(str2)>TH
    ii = strfind(str2,' ');
    ii = ii(find(ii<TH,1,'last'));
    str3 = str2(1:ii);
    str4 = str2(ii+1:end);
    text(0,b-0.08,str3,'units','normalized','FontSize',fs)
    text(0,b-0.16,str4,'units','normalized','FontSize',fs)
    
else
    text(0,b-0.08,str2,'units','normalized','FontSize',fs)
end
end
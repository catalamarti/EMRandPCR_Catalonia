close all
clear all

load('./Data/names.mat')
f = dir('./Figures&Tables/Appendix_EPG_HID/*.png');

for k = 1:8
    
    i1 = 6*(k-1)+1;
    i2 = 6*(k-1)+2;
    i3 = 6*(k-1)+3;
    i4 = 6*(k-1)+4;
    i5 = 6*(k-1)+5;
    i6 = 6*k;
    F1 = imread([f(i1).folder '/' f(i1).name]);
    F2 = imread([f(i2).folder '/' f(i2).name]);
    if k < 8
        F3 = imread([f(i3).folder '/' f(i3).name]);
        F4 = imread([f(i4).folder '/' f(i4).name]);
        F5 = imread([f(i5).folder '/' f(i5).name]);
        F6 = imread([f(i6).folder '/' f(i6).name]);
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
    lx = 0.7;
    rx = 0.7;
    mx = 1;
    fs = 8;
    b = 0.05;
    
    LY = (29.7-uy-2*my-dy)/3;
    LX = (21-lx-rx-mx)/2;
    
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
    
    print(h,['./Figures&Tables/Appendix_EPG_HID/Fig' num2str(k,'%02u') '.pdf'],'-fillpage','-dpdf','-r0')
    
end

delete('./Figures&Tables/Appendix3.pdf')
append_pdfs('./Figures&Tables/Appendix3.pdf','./Figures&Tables/Appendix_EPG_HID/Fig01.pdf',...
    './Figures&Tables/Appendix_EPG_HID/Fig02.pdf','./Figures&Tables/Appendix_EPG_HID/Fig03.pdf',...
    './Figures&Tables/Appendix_EPG_HID/Fig04.pdf','./Figures&Tables/Appendix_EPG_HID/Fig05.pdf',...
    './Figures&Tables/Appendix_EPG_HID/Fig06.pdf','./Figures&Tables/Appendix_EPG_HID/Fig07.pdf',...
    './Figures&Tables/Appendix_EPG_HID/Fig08.pdf')

function [] = caption(i,b,fs,AGA)
TH = 70;
str = ['{\bfFigure ' num2str(i) '}. Daily measures of EMR (maroon) and PCR-based EPGs (black), number of hospitalizations (green), ICU occupancy (red) and mortality (blue) in ' AGA{i} ' over time from late June until the end of August 2020. Hospitalizations, ICU and mortality are averaged over the previous 7 day period.'];
STR = {};
k = 1;
while(length(str)>TH)
    if k == 1
        th = TH;
        k = 0;
    else
        th = TH;
    end
    ii = strfind(str,' ');
    ii = ii(find(ii<th,1,'last'));
    STR = [STR; {str(1:ii)}];
    str = str(ii+1:end);
end
STR = [STR; {str}];
for k = 1:length(STR)
    text(0.03,b-0.05*(k-1),STR{k},'units','normalized','FontSize',fs)
end
end
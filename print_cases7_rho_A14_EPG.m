close all
clear all

STR = {'EPG_HID'};
%'cases7','rho','A14','EPG',
for ii = 1:length(STR)
    
    f = dir(['./Figures&Tables/Appendix_' STR{ii} '/' STR{ii} '_*']);
    
    for k = 1:6
        
        F1 = imread([f(8*(k-1)+1).folder '/' f(8*(k-1)+1).name]);
        F2 = imread([f(8*(k-1)+2).folder '/' f(8*(k-1)+2).name]);
        F3 = imread([f(8*(k-1)+3).folder '/' f(8*(k-1)+3).name]);
        F4 = imread([f(8*(k-1)+4).folder '/' f(8*(k-1)+4).name]);
        if k < 6
            F5 = imread([f(8*(k-1)+5).folder '/' f(8*(k-1)+5).name]);
            F6 = imread([f(8*(k-1)+6).folder '/' f(8*(k-1)+6).name]);
            F7 = imread([f(8*(k-1)+7).folder '/' f(8*(k-1)+7).name]);
            F8 = imread([f(8*(k-1)+8).folder '/' f(8*(k-1)+8).name]);
        end
        
        h=figure(1);
        clf
        K = 1;
        set(h,'Units','centimeters');
        set(h,'Position',[0 0 21/K 29.7/K])
        %set(h,'PaperPositionMode','Auto')
        %set(h,'PaperSize',[21 29.7])
        
        HY = 0.02;
        hy = 0;
        hx = 0;
        
        LY = (1-HY-hy)/4;
        LX = 0.5-hx;
        
        if k<6
            axes('Units','normalized','Position',[hx hy LX LY])
            imagesc(F7)
            axis equal
            axis off
            axes('Units','normalized','Position',[hx+LX hy LX LY])
            imagesc(F8)
            axis equal
            axis off
            axes('Units','normalized','Position',[hx hy+LY LX LY])
            imagesc(F5)
            axis equal
            axis off
            axes('Units','normalized','Position',[hx+LX hy+LY LX LY])
            imagesc(F6)
            axis equal
            axis off
        end
        axes('Units','normalized','Position',[hx hy+LY*2 LX LY])
        imagesc(F3)
        axis equal
        axis off
        axes('Units','normalized','Position',[hx+LX hy+LY*2 LX LY])
        imagesc(F4)
        axis equal
        axis off
        axes('Units','normalized','Position',[hx hy+LY*3 LX LY])
        imagesc(F1)
        axis equal
        axis off
        axes('Units','normalized','Position',[hx+LX hy+LY*3 LX LY])
        imagesc(F2)
        axis equal
        axis off
        
        print(h,['./Figures&Tables/Appendix_' STR{ii} '/Fig' num2str(k,'%02u') '.pdf'], '-fillpage','-dpdf','-r0')
        
    end
    
    append_pdfs(['./Figures&Tables/Appendix_' STR{ii} '.pdf'],['./Figures&Tables/Appendix_' STR{ii} '/Fig01.pdf'],...
        ['./Figures&Tables/Appendix_' STR{ii} '/Fig02.pdf'],['./Figures&Tables/Appendix_' STR{ii} '/Fig03.pdf'],...
        ['./Figures&Tables/Appendix_' STR{ii} '/Fig04.pdf'],['./Figures&Tables/Appendix_' STR{ii} '/Fig05.pdf'],...
        ['./Figures&Tables/Appendix_' STR{ii} '/Fig06.pdf']);
    
end

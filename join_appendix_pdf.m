close all
clear all

STR = {'cases7','rho','A14','EPG','EPG_HID'};

for ii = 1:length(STR)
    
    append_pdfs(['./Figures&Tables/Appendix_' STR{ii} '.pdf'],['./Figures&Tables/Appendix_' STR{ii} '/Fig01.pdf'],...
        ['./Figures&Tables/Appendix_' STR{ii} '/Fig02.pdf'],['./Figures&Tables/Appendix_' STR{ii} '/Fig03.pdf'],...
        ['./Figures&Tables/Appendix_' STR{ii} '/Fig04.pdf'],['./Figures&Tables/Appendix_' STR{ii} '/Fig05.pdf'],...
        ['./Figures&Tables/Appendix_' STR{ii} '/Fig06.pdf']);
    
end
function [ A ] = Bar3Plot2Tikz2Pdf2Png( F,bra,Folder,file,size )

    filename=Bar3Tikz(F,bra,fullfile(Folder,file));
    cd(Folder)
    A=lualatex2pdf2preview( fullfile(pwd,file),0,600 );
    A=imresize(A,[size,NaN]);
    cd('..')
end


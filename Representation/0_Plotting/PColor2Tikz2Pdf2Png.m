function [ A ] = PColor2Tikz2Pdf2Png( bra,ket,indexes,Folder,file,size )
    
    filename=PColorTikz(bra,ket,indexes,fullfile(Folder,file));
    cd(Folder)
    A=lualatex3pdf2preview( fullfile(pwd,file),0,600 );
    A=imresize(A,[size,NaN]);
    cd('..')
end


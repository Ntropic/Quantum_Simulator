function [ A ] = Plot2Tikz2Pdf2Png( Folder,file,size )

    filename=Plot2Tikz(fullfile(Folder,file));
    cd(Folder)
    A=lualatex2pdf2preview( file,0,600 );
    A=imresize(A,[size,NaN]);
    cd('..')
end
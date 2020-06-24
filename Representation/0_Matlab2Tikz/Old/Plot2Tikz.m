function [ filename ] = Plot2Tikz( name )
%BAR3TIKZ creates a tikz file from a barplot in figure,axes F={f,ax}

filename=[name,'.tex'];
matlab2tikz(filename,'mathmode',true,'parseStrings',false);

fid=fopen(filename,'r+');
lines=fgetl(fid);
while lines(1)=='%'
    lines=fgetl(fid);
end
%Get text
text={lines};
while ischar(lines)
        lines=fgetl(fid);
        text={text{:},lines};
end
fclose(fid);
%Pagecolor for beamer format ->black!2
header={'\documentclass[border={25pt 10pt 25pt 10pt}]{standalone}','\usepackage{xcolor}','\pagecolor{black!2}','\usepackage{tikz}','\usepackage{pgfplots}','\usepackage{braket}','\usepgfplotslibrary{patchplots,colormaps}','\usetikzlibrary{matrix,decorations,calc, positioning,arrows,backgrounds,fit,plotmarks}','\begin{document}'};
foot='\end{document}';
while ischar(text{end})==0
    text={text{1:end-1}};
end
text={header{:},text{:},foot};
text=strjoin(text,'\n');
fid=fopen(filename,'w');
fwrite(fid,text);
fclose(fid);
end
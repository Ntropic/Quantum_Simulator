function [ filename ] = PColor2Tikz( name,mode2 )
%BAR3TIKZ creates a tikz file from a barplot in figure,axes F={f,ax}
if nargin==1
    mode2='';
end

filename=[name,'.tex'];
matlab2tikz(filename,'mathmode',true,'parseStrings',false);

fid=fopen(filename,'r+');
lines=fgetl(fid);
while lines(1)=='%'
    lines=fgetl(fid);
end
%Get text
text={lines};
if any(strfind(mode2,'noedge'))==0
    while ischar(lines)
        lines=fgetl(fid);
        if ischar(lines)
            if any(strfind(lines,'title'))
                line2='title style={at={(0.5,-0.05)},anchor=north,yshift=-0.1},';
                text={text{:},line2};
            end
        end
        text={text{:},lines};
    end
else
    while ischar(lines)
        lines=fgetl(fid);
        if ischar(lines)
            if any(strfind(lines,'\addplot['))
                line2='\definecolor{linecolor}{rgb}{0.2081,0.1663,0.5292}';
                text={text{:},line2,lines};
                found=0;
                while ischar(lines) && found==0
                    if ischar(lines)
                        lines=fgetl(fid);
                        if any(strfind(lines,'draw=black'))
                            lines=strrep(lines,'black','linecolor');
                        end
                        text={text{:},lines};
                    end
                end
            end
            if any(strfind(lines,'title'))
                line2='title style={at={(0.5,-0.05)},anchor=north,yshift=-0.1},';
                text={text{:},line2};
            end
        end
        text={text{:},lines};
    end
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
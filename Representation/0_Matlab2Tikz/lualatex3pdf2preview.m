function [ A,h ] = lualatex3pdf2preview( file_name,plotting,res,h )
%TEX2PDF2PREVIEW creates pdf's from .tex files and creates png's of the pdf
%to preview the result. -> executes within folder
if nargin==1
    res=600;
	plotting=1;
    h=figure();
elseif nargin==2
    res=600;
    if plotting==1
        h=figure();
        
    end
elseif nargin==3
    if plotting==1
        h=figure();
    end
end

[filepath,name,ext]=fileparts(file_name);

p=pwd;
if any(strfind(filepath,'\'));
    elem=strsplit(filepath,'\');
else
    elem=strsplit(filepath,'/');
end

%Copy files to name folder

[status,cmdout]=system(['lualatex ',name,'.tex']);
if status~=0 && status~=1
    fprintf(cmdout)
    error(['Conversion to pdf was unsuccessful.'])
end
movefile([name '.pdf'],[name '_re.pdf']);
[status,cmdout]=system(['pdftocairo -png -r ',num2str(res),' ',name,'_re.pdf']);
if status==1
    fprintf(cmdout)
    error(['Conversion from pdf to png was unsuccessful.'])
end

A=imread([name,'_re-1.png']);
if plotting==1
    B=imresize(A,[200,NaN]);
    imshow(B);
    axis equal;
    %ax=get(h,'children');
    %drawnow;
    %pixel_position=getpixelposition(ax);
end
if plotting==0
    h=[];
end
end


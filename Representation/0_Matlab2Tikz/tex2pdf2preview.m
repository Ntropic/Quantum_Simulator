function [ A,h ] = tex2pdf2preview( file_name,plotting,res,h )
%TEX2PDF2PREVIEW creates pdf's from .tex files and creates png's of the pdf
%to preview the result.
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
filepath=fullfile(elem{1:end});
if length(filepath)>0
    cd(filepath);
end

[status,cmdout]=system(['lualatex --aux-directory=auxilliary_files ',name,'.tex']);
if status~=0 && status~=1
    fprintf(cmdout)
    error(['Conversion to pdf was unsuccessful.'])
end

[status,cmdout]=system(['pdftocairo -png -r ',num2str(res),' ',name,'.pdf']);
if status==1
    fprintf(cmdout)
    error(['Conversion from pdf to png was unsuccessful.'])
end

A=imread([name,'-1.png']);
if plotting==1
    A=imresize(A,[200,NaN]);
    imshow(A);
    axis equal;
end
if plotting==0
    h=[];
end

cd(p);
end


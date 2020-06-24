function [ A,h ] = lualatex2pdf2preview( file_name,plotting,res,h )
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
filepath2=fullfile(filepath,name);
if exist(filepath2)==0 || exist(filepath2)==2
    mkdir(name);
end
cd(filepath2);

%Copy files to name folder

[status,cmdout]=system(['lualatex ',fullfile('..',name),'.tex']);
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
%Move pdf and png's to previous directory
d=dir();
for i=3:length(dir)
    if any(findstr(d(i).name,name)) && any(findstr(d(i).name,'.png'))
        movefile(d(i).name,fullfile('..',d(i).name));
    end
end
movefile([name '.pdf'],fullfile('..',[name '.pdf']));
cd(p);
end


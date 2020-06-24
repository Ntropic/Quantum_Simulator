function [ filename, M_pre,M ] = Mat2Tex( M,mode,file_name )
%Mat2TeX creats LaTeX matrices (for quick and dirty matrices use latex(M) 
%        instead) with a prefactor
%extra: (default) maximizes number of ones 
%       'min_size'  - Minimizes size of largest matrix element (slow)
%        'no_opt'   - No optimization
%
%       'vpa'       - to replace fractions by decimals
%       'no_head'   - no header
%       'no_math'   - no math environment
%       'alt_math'  - different math environment (if no header) \[
%       'math'      - default math environment (if no header) \begin{equation}

if nargin==1
    mode='';
end

%% Create Matrix and Prefactor
%Find the ideal prefactor to get out of the matrix 
M_one=length(find(M==1))*ones(size(M));
if length(strfind(mode,'min_size'))==0
    for i=1:size(M,1)
        for j=1:size(M,2)
            if M(i,j)~=0
                M_one(i,j)=length(find(M/M(i,j)==1));
            end
        end
    end
    [a ind]=max(M_one(:));
    [i,j]=ind2sub(size(M),ind);
else
    M_lin=simplifyFraction(simplify(M(:)));
    ind_lin=1:length(M_lin);
    ind_lin=ind_lin(M_lin~=0);
    M_lin=M_lin(ind_lin);
    M_test=M_lin;
    M_length=zeros(length(M_lin),1);
    for i=1:length(M_lin);
        M_test=M_test/M_lin(i);
        
        M_char=char(simplifyFraction(simplify(M_test)));
        M_char=M_char(10:end-3);
        M_char_array=strsplit(M_char,'], [');
        M_length(i)=max(cellfun(@(x) length(x),M_char_array));
    end
    %Find minimal length
    [mini ind]=min(M_length);
    [i,j]=ind2sub(size(M),ind_lin(ind));
end

if length(strfind(mode,'no_opt'))~=0
    M_pre=1;
    M=sym(M);
elseif length(strfind(mode,'vpa'))==0
    M_pre=simplifyFraction(simplify(sym(M(i,j)),10));
    M=simplifyFraction(simplify(sym(M/M_pre),10));
else
    M_pre=vpa(simplifyFraction(simplify(sym(M(i,j))),10));
    M=vpa(simplifyFraction(simplify(sym(M/M_pre)),10));
end

%% Create LaTeX
corrected_latex=latex(M);
cs(1:size(M,2))='c';
corrected_latex=corrected_latex(strfind(corrected_latex,cs)+1+length(cs):end-18);
corrected_latex=['\begin{pmatrix}\n' corrected_latex '\n\end{pmatrix}'];
if M_pre~=1
    b=sym('b');
    prefactor=latex(b*(M_pre));
    prefactor=strrep(prefactor,'\, b\,','');
    prefactor=strrep(prefactor,'\, b','');
    prefactor=strrep(prefactor,'\,b','');
    prefactor=strrep(prefactor,'b\,','');
    tex_code=[prefactor,'\cdot',corrected_latex];
else
    tex_code=[corrected_latex];
end
tex_code=strrep(tex_code,'\\','\\\n');
tex_code=strrep(tex_code,'\','\\');
tex_code=strrep(tex_code,'\\n','\n');

tex_start='\\documentclass[border={25pt 10pt 25pt 10pt}]{standalone} \n\\usepackage{amsmath} \n\n\\begin{document}\n';
tex_end='\\end{document}';

if length(strfind(mode,'no_head'))==0
    if length(strfind(mode,'no_math'))==0
        tex_current=[tex_start '$\n' tex_code '\n$\n' tex_end];
    else
        tex_current=[tex_code];
    end
else %No header
    if length(strfind(mode,'no_math'))==0
    	if length(strfind(mode,'alt_math'))==0
            tex_current=['\\begin{equation}\n' tex_code '\n\\end{equation}\n'];
        else
            tex_current=['\\[\n' tex_code '\n\\]\n'];
        end
    else
        tex_current=[tex_code];
    end
end

%% Save Files
curr_folder=pwd;
save_folder=fullfile(curr_folder,'Matrices');
if exist(save_folder)==0
    mkdir('Matrices');
end


save_folder2=fullfile(save_folder,file_name);
if exist(save_folder2)==0
    mkdir(save_folder2);
end
fileID=fopen([fullfile(save_folder2,file_name),'.tex'],'w');
fprintf(fileID,tex_current);
fclose(fileID);
filename=fullfile(save_folder2,file_name);
end


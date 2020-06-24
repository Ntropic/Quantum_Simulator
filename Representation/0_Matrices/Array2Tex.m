function [ tex_current ] = Array2Tex( array,mode,bra,ket )
%ARRAY2TEX transforms a string array into a latex array
%       'no_head'   - no header
%       'no_environment
%       'none'      - make pmatrix
%       'math'      - not inside math environment
%       'alt_math'  - within a different math environment
%       'tabular'   - inside tabular environment tabular
%       'var_width' - variable width for every column
%       'tikz'      - Create tikz Matrix with brakets (sides) around matrix
%                   -> Add bra,ket for tikz

if nargin==1
    mode=cell(1);
    mode{1}='';
end

tex_start='\\documentclass[border={25pt 10pt 25pt 10pt}]{standalone} \n\\usepackage{rotating} \n\\usepackage{braket} \n\\usepackage{amsmath} \n\\usepackage{amssymb} \n\\usepackage{tikz} \n\\usetikzlibrary{matrix,decorations.pathreplacing, calc, positioning}\n\\usepackage{xcolor}\n\n\\begin{document}\n';
tex_end='\n\\end{document}';

%% Create bulk string (tex_code)
tex_code='';

for i=1:size(array,1)
    for j=1:size(array,2)
        %Assume array is a string array -> other formats might be added
        %later
        if length(strfind(mode{1},'tikz'))==0
            tex_code=[tex_code array{i,j} '&'];
        else
            if length(array{i,j})==0
                tex_code=[tex_code '\\phantom{0} &']; 
            else
                tex_code=[tex_code array{i,j} '&']; 
            end
        end
    end
    tex_code(end)=[];
    if i<(size(array,2))
        tex_code=[tex_code '\\\\\n'];
    else
        if length(strfind(mode{1},'tikz'))==0
            tex_code=[tex_code '\n'];
        else
            tex_code=[tex_code '\\\\\n'];
    end
end

%% Add head, environment etc
c_rep(1)='r';
c_rep(2:size(array,2))='c';
if length(strfind(mode{1},'no_head'))==0
    if length(strfind(mode{1},'none'))==1
        tex_current=[tex_start '$\\begin{pmatrix}\n',tex_code,'\n\\end{pmatrix}$' tex_end];
    elseif length(strfind(mode{1},'tabular'))==1
        tex_current=[tex_start '\\begin{tabular}{',c_rep,'}\n',tex_code,'\n\\end{tabular}' tex_end];
    elseif length(strfind(mode{1},'tikz'))>0
        tex_current=[tex_start '\n\\begin{tikzpicture}\n\\matrix [matrix of math nodes,left delimiter=(,right delimiter=)] (m) {\n' tex_code '};\n']; %,row sep=0.5cm,column sep=0.5cm
        %Add nodes
        if nargin>=3
            for i=1:length(bra)
                tex_current=[tex_current '\\node[left=12pt of m-' num2str(i) '-1] (left-' num2str(i) ') {$' bra{i} '$};\n'];
            end
        end
        if nargin==4
            for i=1:length(ket)
                tex_current=[tex_current '\\node[above=5pt of m-1-' num2str(i) '] (top-' num2str(i) ') {$' ket{i} '$};\n'];
            end
        end
        tex_current=[tex_current '\n\\end{tikzpicture}' tex_end];
    elseif length(strfind(mode{1},'no_math'))==0
        tex_current=[tex_start '$\n\\begin{array}{' c_rep '}\n' tex_code '\\end{array}\n$\n' tex_end];
    else
        tex_current=[tex_start tex_code tex_end];
    end
else %No header
    if length(strfind(mode{1},'none'))==1
        if length(strfind(mode{1},'alt_math'))==0
            tex_current=['\\begin{equation}\\begin{pmatrix}\n',tex_code,'\n\\end{pmatrix}\\end{equation}'];
        else
            tex_current=['\\[\\begin{pmatrix}\n',tex_code,'\n\\end{pmatrix}\\]'];
        end
    elseif length(strfind(mode{1},'no_math'))==0
        if length(strfind(mode{1},'tikz'))>0
            tex_current=['\n\\begin{tikzpicture}\n\\matrix [matrix of math nodes,left delimiter=(,right delimiter=)] (m) {\n' tex_code '};\n']; %,row sep=0.5cm,column sep=0.5cm
            %Add nodes
            if nargin>=3
                for i=1:length(bra)
                    tex_current=[tex_current '\\node[left=12pt of m-' num2str(i) '-1] (left-' num2str(i) ') {$' bra{i} '$};\n'];
                end
            end
            if nargin==4
                for i=1:length(ket)
                    tex_current=[tex_current '\\node[above=5pt of m-1-' num2str(i) '] (top-' num2str(i) ') {$' ket{i} '$};\n'];
                end
            end
            tex_current=[tex_current '\n\\end{tikzpicture}' tex_end];
        elseif length(strfind(mode{1},'alt_math'))==0
            tex_current=['\\begin{equation}\n\\begin{array}{' c_rep '}\n' tex_code '\\end{array}\n\\end{equation}\n'];
        else
            tex_current=['\\[\n\\begin{array}{' c_rep '}\n' tex_code '\\end{array}\n\n\\]\n'];
        end
    else
        tex_current=[tex_code];
    end
end
end


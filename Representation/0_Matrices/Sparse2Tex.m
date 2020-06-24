function [ filename ] = Sparse2Tex( S,file_name,mode,varargin)
%SPARSE2TEX creates a square matrix (M) from a sparse matrix (S) that
%contains only the subspace of S and expm(S) 
%Modes: 'none'      - no additionsv -> pmatrix (round braces)
%       'index'     - add indices
%       'fock'      - add fock number (or fock_bin for binary fock numbers)
%       'fock_dense' - Using vertical bra vectors
%Additional:
%       'unit'          - S as a unitary matrix (only non trivial elements)
%       'no_zeros'      - remove zeros
%       'vpa'           - Use doubles instead of fractions
%       'full'          - print full matrix (default)
%       'def_sub'       - use default substitutions for unitary matrices
%
%     	'var_width'     - variable width for every column
%       Future: 'fixed_width'   - same width for all columns
%       'no_head'       - no header
%       'no_math'       - no math environment
%
%mode is either char or a cell (if mode is 'fock' it needs to be a cell
%with mode{2}=n (maximum number of photons per mode) and mode{3}=N (number
%of modes)
%ancilla qbits can be added via mode{4}
%
%index_list,substitution_list,substitution_compare and def_sub are optional
%input arguments
%       index_list           - specifies the indexes that should be printed
%       substitution_list    - specifies the substitute chars (tex) for the
%                              sequences, that are being compared with
%       substitution_compare - if sequence is found it gets substituted by 
%                              the corresponding substitution_list element 

tex_start='\\documentclass[border={25pt 10pt 25pt 10pt}]{standalone} \n\\ usepackage{rotating} \n\\usepackage{braket} \n\\usepackage{amsmath} \n\n\\begin{document}\n';
tex_end='\\end{document}';

eps=1e-6;
%% Prepare data (Matrix entries)
if nargin==1
    mode='';
end

if iscell(mode)==0
    %transform to cell
    mode2=mode;
    mode=cell(1);
    mode{1}=mode2;
end
if strfind(mode{1},'fock')
    if length(mode)<3 && length(mode)>1
        error('Mode must also specify the maximum number of particles per mode (n) and the number of modes N.')
    elseif length(mode)==1
        mode{2}=1;
        mode{3}=log2(size(S,1));
    end
end
do_index_list=0;
do_substitution_list=0;
do_substitution_compare=0;
if  nargin>3
    if mod(nargin,2)~=1
        error('Additional input arguments must be specified via a specification string and the variables.')
    end
    for i=1:(nargin-3)/2
        j=2*i-1;
        k=2*i;
        if any(strfind(varargin{j},'index_list'))==1
            do_index_list=1;
            index_list=varargin{k};
        end
        if any(strfind(varargin{j},'substitution_list'))==1
            do_substitution_list=1;
            substitution_list=varargin{k};
        end   
        if any(strfind(varargin{j},'substitution_compare'))==1
            do_substitution_compare=1;
            substitution_compare=varargin{k};
        end   
    end
end
if do_index_list==0
    if length(strfind(mode{1},'unit'))==1
        %Remove numeric uncertainty
        a=find(S);
        b=find(abs(S)<10^-15);
        [~,pos]=intersect(a,b);
        S(a(pos))=0;
        %Find non trivial values
        [i,j,t]=find(S);
        k=i-j;
        k=find(k);  %Indexes that are not on the diagonal
        ko=[i(k);j(k)];
        
        diag_find=find(i-j==0);
        S_diag=S(i(diag_find)+size(S,1)*(i(diag_find)-1));
        is_one=[];
        if isa(S_diag,'sym')
            for l=1:length(S_diag)
                if length(symvar(S_diag(l)))>0
                    is_one=[is_one,l];
                elseif abs(S_diag(l)-1)>10^-13
                    is_one=[is_one,l];
                end
            end
        else %Double
            for l=1:length(S_diag)
                if abs(S_diag(l)-1)>10^-13
                    is_one=[is_one,l];
                end
            end
        end
        k=[ko(:);is_one'];
        k=sort(unique(k));
    elseif length(strfind(mode{1},'full'))==1
        i=1:length(S);
        k=[i];
        k=sort(unique(k));
    else
        [i,j,t]=find(S);
        k=[i;j];
        k=sort(unique(k));
    end
else
    i=index_list;
    j=index_list;
    k=[i;j];
    k=sort(unique(k));
end
if length(k)==0
    k=1;
end
%Concatenate 
ket_str=k;


%% Prepare Default Substitution List
if length(strfind(mode{1},'def_subs'))~=0
    alpha=sym('alpha');
    beta=sym('beta');
    delta=sym('delta');
    theta=sym('theta');
    assume(alpha,'real');
    assume(beta,'real');
    assume(delta,'real');
    assume(theta,'real');
    subs_comp=[(exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2,- (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 + (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,
               (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2)*1i)/2 - (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2)*1i)/2,(exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp(-(theta*1i)/2))/2 + (exp((alpha*1i)/2)*exp((beta*1i)/2)*exp(delta*1i)*exp((theta*1i)/2))/2,
               (exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp((delta*1i)/2)*exp(-(theta*1i)/2))/2 + (exp(-(alpha*1i)/2)*exp(-(beta*1i)/2)*exp((delta*1i)/2)*exp((theta*1i)/2))/2,- (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp((delta*1i)/2)*exp(-(theta*1i)/2)*1i)/2 + (exp(-(alpha*1i)/2)*exp((beta*1i)/2)*exp((delta*1i)/2)*exp((theta*1i)/2)*1i)/2,
               (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp((delta*1i)/2)*exp(-(theta*1i)/2)*1i)/2 - (exp((alpha*1i)/2)*exp(-(beta*1i)/2)*exp((delta*1i)/2)*exp((theta*1i)/2)*1i)/2,(exp((alpha*1i)/2)*exp((beta*1i)/2)*exp((delta*1i)/2)*exp(-(theta*1i)/2))/2 + (exp((alpha*1i)/2)*exp((beta*1i)/2)*exp((delta*1i)/2)*exp((theta*1i)/2))/2,
               cos(theta/2)*exp(-(alpha*1i)/2-(beta*1i)/2+delta*1i) , sin(theta/2)*exp((alpha*1i)/2-(beta*1i)/2+delta*1i) , 
               -sin(theta/2)*exp(-(alpha*1i)/2+(beta*1i)/2+delta*1i) , cos(theta/2)*exp((alpha*1i)/2 + (beta*1i)/2 + delta*1i) ];

           subs_list={'\hat{U}_{11}','\hat{U}_{12}',
               '\hat{U}_{21}','\hat{U}_{22}',
               '\hat{U}_{11}','\hat{U}_{12}',
               '\hat{U}_{21}','\hat{U}_{22}',
               '\hat{U}_{11}','\hat{U}_{12}',
               '\hat{U}_{21}','\hat{U}_{22}'};
    if do_substitution_list==1
        substitution_list={substitution_list{:};subs_list{:}}';
    else
        substitution_list={subs_list{:}}';
    end
    if do_substitution_compare==1
        substitution_compare={substitution_compare{:};subs_comp(:)};
    else
        substitution_compare=subs_comp(:);
    end
    do_substitution_list=1;
    do_substitution_compare=1;
end

%% Create bra and ket elements
M=sym(full(S(k,k)));

if length(strfind(mode{1},'index'))>=1
    max_length=1;
    bra_str=cell(length(k),1);
    ket_str=cell(length(k),1);
    long_str=num2str(max(k));
    type=['%0' num2str(length(long_str)) 'd'];
    for i=1:length(k)
        ket_str{i}=num2str(k(i),type);
    end
    ket_str=strcat({'\\mathbf{['},ket_str,{']}'});
    bra_str=ket_str;
elseif length(strfind(mode{1},'fock'))>=1
    %Create Fock-states
    max_length=1;
    
    bra_str=bra_fock_tex(k,mode);
    ket_str=ket_fock_tex(k,mode);
    ket_str=strcat({'\\mathbf{'},ket_str,{'}'});
    bra_str=strcat({'\\mathbf{'},bra_str,{'}'});
    if length(strfind(mode{1},'fock_dense'))==1
        [vert_str]=vertfock_tex(k,mode,'mathbf');
        vert_str=strcat({'\\mathbf{'},vert_str,{'}'});
        ket_str=vert_str;
    end
else %if length(strfind(mode{1},'none'))==1
    max_length=1;
    M_mode=0;
    bra_str=cell(length(k),1);
    ket_str=cell(length(k),1);
    for i=1:length(k)
        ket_str{i}='';
        bra_str{i}='';
    end
end

%% Create strings for entries 
%No substitutions
if do_substitution_list==0 || do_substitution_compare==0
    M_str=cell(length(k),length(k));
    if length(strfind(mode{1},'vpa'))==0
        for o=1:length(k)
            for p=1:length(k)
                M_str{o,p}=latex(simplifyFraction(simplify(M(o,p))));
                M_str{o,p}=strrep(M_str{o,p},'\','\\');
            end
        end
    else
        for o=1:length(k)
            for p=1:length(k)
                M_str{o,p}=latex(vpa(simplify(M(o,p))));
                M_str{o,p}=strrep(M_str{o,p},'\','\\');
            end
        end
    end
%Substitutions
elseif do_substitution_compare==1 && do_substitution_list==1 %both needed
    if size(substitution_compare)~=size(substitution_list)
        error('Substitution list and substitution compare list have to be of the same size');
    end
    substitution_compare=[substitution_compare(:)];
    substitution_list={substitution_list{:}};
    
    %Create syms for presubstitution in whole matrix
    sym1=[];
    sym1_tex={};
    sym_length=length(num2str(length(substitution_list)));
    for ind=1:length(substitution_list)
        sym1=[sym1,sym(['sub' num2str(ind,['%0' num2str(sym_length) 'd'])])];
        sym1_tex={sym1_tex{:},latex(sym1(ind))};
    end
    
    %Create array elements
    lens=zeros(length(substitution_list),1);
    for o=1:length(k)
        for p=1:length(k)
            %Try all substitutions
            if length(char(M(o,p)))>length(sym1_tex{1})
                for q=1:length(substitution_list)
                    M(o,p)=subs(M(o,p),substitution_compare(q),sym1(q));
                    test=char(M(o,p));
                    lens(q)=length(test);
                end
                [ala minind]=min(lens);
                M(o,p)=subs(M(o,p),substitution_compare(minind),sym1(minind));
                M(o,p)=subs(M(o,p),substitution_compare,sym1');
            end
            if length(strfind(mode{1},'vpa'))==0
                M_str{o,p}=latex(simplifyFraction(simplify(M(o,p))));
            else
                M_str{o,p}=latex(vpa(simplify(M(o,p))));
            end
            M_str{o,p}=strrep(M_str{o,p},'\','\\');
        end
    end
    
    %Replace substitutes by corresponding substitution_list elements
    for ind=1:length(substitution_list)
        M_str=cellfun(@(x) strrep(x,sym1_tex(ind),substitution_list(ind)),M_str);
    end
    
    %Replace zeros with lightgray zeros
    M_str=cellfun(@(x) strrep(x,{'0'},{'\\textcolor{lightgray}{0}'}),M_str);
end
if length(strfind(mode{1},'no_zeros'))>0
    for o=1:length(k)
        for p=1:length(k)
            M_str{o,p}=strrep(M_str{o,p},'0','');
        end
    end    
end

%% Create Array
%Top Line 
if length(strfind(mode{1},'none'))==0
    bra_str={'',bra_str{:}}';
    array=[bra_str,[ket_str';M_str]];
else
    array=M_str;
end

%% Create string
%% Create string
if length(strfind(mode{1},'tikz'))==0
    if length(strfind(mode{1},'none'))==0
        bra_str={'',bra_str{:}}';
        array=[bra_str,[ket_str';M_str]];
    else
        array=M_str2;
    end
    tex_current=Array2Tex(array,mode);
else
    tex_current=Array2Tex(M_str,mode,{bra_str{2:end}},ket_str);
end

%% Save Files
curr_folder=pwd;
save_folder=fullfile(curr_folder,'Sparse2Square_Matrices');
if exist(save_folder)==0
    mkdir('Sparse2Square_Matrices');
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
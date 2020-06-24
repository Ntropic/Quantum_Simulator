function [ M ket_str M_mode ] = Sparse2Square( S,mode,varargin)
%SPARSE2SQUARE creates a square matrix (M) from a sparse matrix (S) that
%contains only the subspace of S and expm(S)
%Modes: 'none'  - no additions
%       'index' - add indices
%       'fock'  - add fock number (or fock_bin for binary fock numbers)
%       'fock_dense' - Using vertical bra vectors
%       'fock_bin'   - 
%Additional:
%       'unit'      - Treat S as a unitary matrix
%       'full'      - print full matrix
%       'free_line' - adds a free line between top row and matrix
%       'def_sub'   - use default substitutions for unitary matrices
%     One Day?  'var_width' - variable width of the rows
%
%mode is either char or a cell (if mode is 'fock' it needs to be a cell
%with mode{2}=n (maximum number of photons per mode) and mode{3}=N (number
%of modes)
%ancilla qbits can be added via mode{4}
%
%index_list,substitution_list,substitution_compare and def_sub are optional
%input arguments
%       index_list           - specifies the indexes that should be printed
%       substitution_list    - specifies the substitute chars for the
%                              sequences, that are being compared with
%       substitution_compare - if sequence is found it gets substituted by 
%                              the corresponding substitution_list element 
eps=1e-6;
%Prepare data
if nargin==1
    mode='none';
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
if  nargin>2
    if mod(nargin,2)~=0
        error('Additional input arguments must be specified via a specification string and the variables.')
    end
    for i=1:nargin/2-1
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
        [i,j,t]=find(S); %Only use dimensions with off diagonal elements
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
    subs_list={'U11','U12',
               'U21','U22',
               'U11','U12',
               'U21','U22',
               'U11','U12',
               'U21','U22'};
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

%% Create new matrix 
M=full(S(k,k));

max_length=1;
if any(strfind(mode{1},'none'))>1
    max_length=1;
    M_mode=0;
    bra_str=cell(length(k),1);
    ket_str=cell(length(k),1);
    for i=1:length(k)
        ket_str{i}='';
        bra_str{i}='';
    end
elseif any(strfind(mode{1},'index'))>1
    max_length=1;
    bra_str=cell(length(k),1);
    ket_str=cell(length(k),1);
    long_str=num2str(max(k));
    type=['%0' num2str(length(long_str)) 'd'];
    for i=1:length(k)
        ket_str{i}=num2str(k(i),type);
    end
    ket_str=strcat({'['},ket_str,{']'});
    bra_str=ket_str;
elseif any(strfind(mode{1},'fock'))>0
    %Create Fock-states
    max_length=1;
    bra_str=bra_fock_str(k,mode);
    ket_str=ket_fock_str(k,mode);
    if length(strfind(mode{1},'fock_dense'))==1
        [vert_str,~,~,fock_comp_len]=vertfock_str(k,mode);
    end
end

%Create strings for entries 
if isa(M,'sym')==0 %Not symbolic
    M_str=cell(length(k));
    if max(abs(imag(M(:))))>eps
        for o=1:length(k)
            for p=1:length(k)
                if abs(M(o,p))>eps
                    real_phi=sprintf('%0.2f',full(real(M(o,p))));
                    imag_phi=sprintf('%0.2f',full(imag(M(o,p))));
                    if real_phi(1)~='-'
                        real_phi=[' ' real_phi];
                    end
                    if imag_phi(1)~='-'
                        imag_phi=['+i' imag_phi];
                    else
                        imag_phi=['-i' imag_phi(2:end)];
                    end
                    M_str{o,p}=[real_phi imag_phi];
                    if max_length<length(M_str{o,p})
                        max_length=length(M_str{o,p});
                    end
                else
                    M_str{o,p}='0';
                end
            end
        end
    else
        %Only real part
        for o=1:length(k)
            for p=1:length(k)
                real_phi=sprintf('%0.2f',full(real(M(o,p))));
                if abs(full(real(M(o,p)))-round(full(real(M(o,p)))))<eps
                    real_phi=sprintf('%d',full(round(real(M(o,p)))));
                end
                if real_phi(1)~='-'
                    real_phi=[' ' real_phi];
                end
                M_str{o,p}=[real_phi];
                if max_length<length(M_str{o,p})
                    max_length=length(M_str{o,p});
                end
            end
        end
    end
else %is symbolic Use: [ string ] = sym2str( symbolic_matrix )
    if do_substitution_list==0
        M_str=cell(length(k));
        for o=1:length(k)
            for p=1:length(k)
                M_str{o,p}=char(simplifyFraction(simplify(M(o,p))));
                if max_length<length(M_str{o,p})
                   max_length=length(M_str{o,p});
                end
            end
        end
    elseif do_substitution_compare==0
        M_str=cell(length(k));
        M_prior=char(M);
        M_prior=M_prior(10:end-2);
        M_vec=strsplit(M_prior,'], [');
        M_str=regexp(M_vec,', ','split');
        M_str=vertcat(M_str{:});
        if do_index_list==1
            for o=index_list
                for p=index_list
                    M_str{o,p}=substitution_list{mod(find(k==o)+1,length(substitution_list))+1,mod(find(k==p)+1,length(substitution_list))+1};
                end
            end
        end
        M_str=M_str(k,k);
        M_str_len=max(max(cellfun(@length,M_str)));
        if max_length<M_str_len
            max_length=M_str_len;
        end
    elseif do_substitution_compare==1
        if size(substitution_compare)~=size(substitution_list)
            error('Substitution list and substitution compare list have to be of the same size');
        end
        M_str=cell(length(k));
        
        M_prior=char(M);
        M_prior=M_prior(10:end-3);
        M_vec=strsplit(M_prior,'], [');
        M_str=regexp(M_vec,', ','split');
        M_str=vertcat(M_str{:});
        [ind1 ind2]=find(strcmp(M_str,'1')+strcmp(M_str,'0')==0);
        
        subs_cmp_str=cell(size(substitution_compare));
        for p=1:size(substitution_compare,1)
            for q=1:size(substitution_compare,2)
                subs_cmp_str{p,q}=char(substitution_compare(p,q));
            end
        end
        
        for ind=1:length(ind1)
            o=ind1(ind);
            p=ind2(ind);
            found=0;
            for q=1:size(substitution_compare,1)
                for r=1:size(substitution_compare,2)
                    if strcmp(subs_cmp_str{q,r},M_str{o,p})
                        M_str{o,p}=substitution_list{q,r};
                        found=1;
                    end
                end
            end
            if found==0
                for q=1:size(substitution_compare,1)
                    for r=1:size(substitution_compare,2)
                        if isAlways(M(o,p)==substitution_compare(q,r),'Unknown','false')
                            M_str{o,p}=substitution_list{q,r};
                        end
                    end
                end
            end
        end
        M_str_len=max(max(cellfun(@length,M_str)));
        if max_length<M_str_len
            max_length=M_str_len;
        end        
    end
end


if length(strfind(mode{1},'fock_dense'))>0
    lengths_of_ket_index=cellfun('length',ket_str);
    length_of_index_top=fock_comp_len;
else
    lengths_of_ket_index=max(cellfun('length',ket_str));
    length_of_index_top=max(cellfun('length',bra_str));
end    
if max_length<length_of_index_top
    max_length=length_of_index_top;
end 
max_length=max_length+2;

    
%Create M_mode
M_mode='';
%First row(s) (multiple if dense parameter was found)
if length(strfind(mode{1},'none'))==0
    if length(strfind(mode{1},'index'))>0   
        len=length(M_mode)+lengths_of_ket_index+2;
        M_mode(len+max_length-length(bra_str{1}):len+max_length+8)=[' <strong>' bra_str{1}];
        len=length(M_mode);        
        for o=2:length(k)
            M_mode(len+max_length-length(bra_str{o})+1:len+max_length)=bra_str{o};
            len=length(M_mode);
        end
        M_mode=[M_mode '</strong>\n'];
    elseif length(strfind(mode{1},'fock_dense'))==0
        len=length(M_mode)+lengths_of_ket_index+2;
        M_mode(len+max_length-length(bra_str{1}):len+max_length+7)=['<strong>' bra_str{1}];
        len=length(M_mode);        
        for o=2:length(k)
            M_mode(len+max_length-length(bra_str{o})+1:len+max_length)=bra_str{o};
            len=length(M_mode);
        end
        M_mode=[M_mode '</strong>\n'];
    else %vertical ket 
        verter=regexp(vert_str,'\\n','split');
        vert_str_split=vertcat(verter{:});
        for i=1:size(vert_str_split,2)
            len=length(M_mode)+lengths_of_ket_index+2;
            M_mode(len+max_length-length_of_index_top:len+max_length+7-length_of_index_top+length(vert_str_split{1,i}))=['<strong>' vert_str_split{1,i}];
            len=length(M_mode);        
            for o=2:length(k)
                M_mode(len+max_length-length_of_index_top+1:len+max_length-length_of_index_top+length(vert_str_split{o,i}))=vert_str_split{o,i};
                len=length(M_mode);
            end
            M_mode=[M_mode '</strong>\n'];
        end
    end
end
if length(strfind(mode{1},'free_line'))>0
    for i=1:length(strfind(mode{1},'free_line'))
        M_mode=[M_mode '\n'];
    end
end

for o=1:length(k)
    M_mode=[M_mode '<strong>'];
    len=length(M_mode);
    M_mode(len+1:len+length(ket_str{o}))=ket_str{o};
    M_mode=[M_mode '</strong>'];
    len=length(M_mode)+lengths_of_ket_index+1-length(ket_str{o});
    for p=1:length(k)
        M_mode(len+max_length-length(M_str{o,p})+1:len+max_length)=M_str{o,p};
        len=length(M_mode);
    end
    M_mode=[M_mode '\n'];
    len=length(M_mode);
end
end
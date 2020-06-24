function [ fun_mat, fun_vars_new ] = fun_mat_subs( fun_mat,fun_vars,fun_var_subs )
%FUN_MAT_SUBS substitutes variables fun_vars with other variables encoded
%in fun_var_subs
% Input:
%   fun_mat         -   is a function handle of a reshaped matrix
%   fun_vars        -   is a cell of variables in fun_mat
%   fun_var_subs    -   is a cell or list of symbolic expressions and 
%                       numeric expressions for the substitution
% Output:
%   fun_mat         -   new function handle of the matrix after
%                       substitution
%   fun_vars        -   list of current variables in the changed function
%                       handle (also of course this list can be found in
%                       the fun_mat, but this is more convenient for
%                       extraction
if length(fun_vars)>0
    if isa(fun_mat,'function_handle')
        if iscell(fun_var_subs)==0
            %Determine remaining parameters - are all symbolic,numeric, or are
            %there both types?
            fun_vars_subs=[];
            c=0;
            d=zeros(1,length(fun_var_subs));
            for i=1:length(fun_var_subs)
                symmian=symvar(fun_var_subs(i));   
                fun_vars_subs=[fun_vars_subs,symmian];
                if length(symmian)>0
                    c=c+1;
                    d(i)=1;
                end
            end

            if c==length(fun_var_subs) %All symbolic expressions (or the subs)
                fun_vars_new=strsplit(char(symvar(fun_vars_subs)),{'(',')','[',']',', '});
                if length(fun_vars_subs)>1
                    fun_vars_new=fun_vars_new(2:end-1);
                end

                [c_mat,c_vars]=fun_mat2str(fun_mat);

                %Change c_vars
                c=strsplit(c_vars,{'(',')'});
                c_vars=[c{1} '(' strjoin(fun_vars_new,',') ')' c{end-1} '('];

                for i=1:length(fun_var_subs)
                    c_mat=strrep(c_mat,fun_vars{i},char(fun_var_subs(i)));
                end
                fun_mat=str2fun_mat(c_mat,c_vars);

            elseif c==0 %All numeric expressions
                b=num2cell(double(fun_var_subs));
                fun_mat=fun_mat(b{:});
                fun_vars_new={};

            else %Both numeric and symbolic expressions -> needs symbolic substitution (avoid at all cost!)
                %First substitute numeric values
                sym_mat=sym(fun_mat);
                for i=1:length(fun_var_subs)
                    if d(i)==0    %Is numeric
                        sym_mat=subs(sym_mat,fun_vars{i},fun_var_subs(i));
                    end
                end
                fun_mat=matlabFunction(sym_mat);

                fun_vars_new=strsplit(char(symvar(fun_vars_subs)),{'(',')','[',']',', '});

                [c_mat,c_vars]=fun_mat2str(fun_mat);

                %Change c_vars
                c=strsplit(c_vars,{'(',')'});
                c_vars=[c{1} '(' strjoin(fun_vars_new,',') ')' c{end-1} '('];

                c=1;
                for i=1:length(fun_var_subs)
                    if d(i)==1
                        c_mat=strrep(c_mat,fun_vars{i},char(fun_var_subs(1)));
                        c=c+1;
                    end
                end
                fun_mat=str2fun_mat(c_mat,c_vars);
            end


        else %Is a cell
            % First: Are all "fun_var_subs(2*k)" parameters symbolic, numeric, or 
            %        are both types present
            % Second: Are all parameter descriptors "fun_var_subs(2*k-1)" - the
            %         parameters to be replaced - simple parameters or functions
            %         themselves? -> Avoid at any cost!!! This is extremely slow!

            if mod(length(fun_var_subs),2)~=0
                error('Cell substitutions need to have length divisible by 2.');
         	end
            
            k=1:length(fun_var_subs)/2;
            fun_var_desc=fun_var_subs(2*k-1);
            fun_var_subs=fun_var_subs(2*k);
            
            
            %First determine types of subs
            c=0;
            d=zeros(1,length(fun_var_subs));
            fun_vars_subs=[];
            for i=1:length(fun_var_subs)
                symmian=symvar(fun_var_subs{i});
                fun_vars_subs=[fun_vars_subs,symmian];
                if length(symmian)>0
                    c=c+1;
                    d(i)=1;
                end
            end
            %Second determine if these are simple substitutions or not!
            simple=zeros(1,length(fun_var_subs));
            for i=1:length(fun_var_subs)
                a=char(fun_var_desc{i});
                b=char(symvar(fun_var_desc{i}));
                simple(i)=strcmp(a,b);
            end

            %Do the substitutions case dependent (3 cases: 
            %   1. all_num simple, 
            %   2. all_sym simple,
            %   3. non simple or (some_num) -> both require 

            if all(simple) && ( c==0 || c==length(fun_var_subs)) %Cases 1. and 2.
                %Reorder parameters
                ind=zeros(1,length(fun_var_desc));
                for i=1:length(fun_var_desc)
                    name=char(fun_var_desc{i});
                    a=find(strcmp(fun_vars,name));
                    if length(a)==1
                        ind(i)=a;
                    end
                end
                findi=find(ind);
                ind=ind(findi);
                [m n]=sort(ind);
                findi=findi(n);
                %First
                if c==0 %All numeric expressions
                    %how_to_order
                    fun_new=zeros(length(ind),1);
                    for p=1:length(ind)
                        fun_new(p)=double(fun_var_subs{findi(p)});
                    end
                    a=num2cell(fun_new);
                    fun_mat=fun_mat(a{:});
                    fun_vars_new={};

                %Second
                elseif c==length(fun_var_subs) %All symbolic expressions
                    fun_vars_new=strsplit(char(symvar(fun_vars_subs)),{'(',')','[',']',', '});
                    if length(fun_vars_subs)>1
                        fun_vars_new=fun_vars_new(2:end-1);
                    end

                    [c_mat,c_vars]=fun_mat2str(fun_mat);

                    %Change c_vars
                    c=strsplit(c_vars,{'(',')'});
                    c_vars=[c{1} '(' strjoin(fun_vars_new,',') ')' c{end-1} '('];

                    for i=1:length(fun_var_subs)
                        if length(symvar(fun_var_subs{i}))>1
                            c_mat=strrep(c_mat,fun_vars{i},['(' char(fun_var_subs{i}) ')']);
                        else
                            c_mat=strrep(c_mat,fun_vars{i},char(fun_var_subs{i}));
                        end
                    end
                    fun_mat=str2fun_mat(c_mat,c_vars);

                end
            else %Either both symbolic and numeric or non simple! (takes long time)
                sym_mat=sym(fun_mat);
                for k=1:length(fun_var_subs)
                    sym_mat=subs(sym_mat,fun_var_desc{i},fun_var_subs(i));
                end

                fun_mat=matlabFunction(sym_mat); 
                [c_mat,c_vars]=fun_mat2str(fun_mat);

                fun_vars_new=strsplit(char(c_vars),{'(',')',', '});
                fun_vars_new=fun_vars_new(2:end-2);
            end
        end
    elseif isa(fun_mat,'sym')
        error('fun_mat_subs only supports function_handles as input.');
    else %Double
        fun_vars_new={};
    end
else %Differentiate cases of different data types
    if isa(fun_mat,'function_handle')
        fun_vars_new=fun_vars;
    elseif isa(fun_mat,'sym')
        error('fun_mat_subs only supports function_handles as input.');
    else %Double
        fun_vars_new={};
    end
end
end


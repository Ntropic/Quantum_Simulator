function [matrix anc_matrix sizes anc_sizes] = Expansion2Matrix(gate,current_param,plotting)
%Expansion2Matrix creates a matrix out of the elementary gate expansion
%defined in exp_steps of the "gate"
%Check if gates have been prepared?!
%       ->  Prepare Gates via "Gates_Tables_Prepare"
if nargin==1
    current_param={};
    plotting=1;
elseif nargin==2
    plotting=1;
end
fac=500;

sizes=gate.size;
anc_sizes=gate.anc_size;
N=sizes+anc_sizes;


index=[gate.exp_steps.index];
matrix=[gate.exp_steps.matrix];
par=[gate.exp_steps.param];
param=unique([gate.exp_steps.param{:}]);
global_phase=gate.exp_steps.global_phase;



if length(current_param)==0  
    if length(param)>0
        warning('Not all parameters have been substituted')
    end
    %Embed Gates iteratively
    if isa(matrix{1},'function_handle')==0
        if length(global_phase)>0
            anc_matrix=Embed_Gate(matrix{1},index{1},N)*global_phase;
        else
            anc_matrix=Embed_Gate(matrix{1},index{1},N);
        end
    else
        if length(global_phase)>0
            anc_matrix=Embed_Gate(sym(matrix{1}),index{1},N)*global_phase;
        else
            anc_matrix=Embed_Gate(sym(matrix{1}),index{1},N);
        end
    end
    nums=floor(length(matrix)/fac);
    lens=length(matrix);
    if nums<1
        for i=2:length(matrix)
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N)*anc_matrix;
            end
        end
    else
        timer=tic();
        for i=2:fac
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N)*anc_matrix;
            end
        end
        t=toc(timer);
        if plotting==1
            fprintf(['    -> Remaining time: ' num2str(lens/fac*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac) '/' num2str(lens) '\n']);
        end
        for j=2:nums
            for i=fac*(j-1)+1:fac*j
                if isa(matrix{i},'function_handle')==0
                    anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
                else
                    anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N)*anc_matrix;
                end
            end
            t=toc(timer);
            if plotting==1
                fprintf(['    -> Remaining time: ' num2str(lens/(fac*j)*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac*j) '/' num2str(lens) '\n']);
            end
        end
        for i=fac*nums+1:length(matrix)
            if isa(matrix{i},'function_handle')==0
                anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
            else
                anc_matrix=Embed_Gate(sym(matrix{i}),index{i},N)*anc_matrix;
            end
        end
        t=toc(timer);
        if plotting==1
            fprintf(['    -> Remaining time: 0s, Elapsed: ' num2str(t) 's, ' num2str(lens) '/' num2str(lens) '\n']);
        end
    end
else %Do last substitutions
    if iscell(current_param)==0   %Is not a cell -> Use variables as input   
        if length(param)~=length(current_param) 
            error('Not all parameters have been substituted. Use Gate2Matrix if you require a symbolic expansion.')
        else
            %Substitute for every elementary gate in list
            for i=1:length(gate.exp_steps.index)
                matrix{i}=fun_mat_subs(matrix{i},par{i},current_param);
            end
        end
    else %Is a cell
        if mod(length(current_param),2)
            error(['Not the right amount of parameter inputs. Has to be divisible by 2 for cell inputs.'])
        else
            %Substitute for every elementary gate in list
            matrix{i}=fun_mat_subs(matrix{i},par{i},current_param);
        end
    end
    %Embed Gates iteratively

    %% Global phase
    %Substitute global phase aswell!!! 
    if length(global_phase)>0
        if iscell(current_param)==0 %Is not a cell -> Use variables as input
            if isa(global_phase,'sym')==1
                pa=strsplit(char(symvar(global_phase)),{'[',']'});
                if length(pa)>1
                    pa=pa{2};
                    pa=strsplit(pa,', ');
                end
            end
            len_pa=length(pa); %Current global_phase parameters
            if len_pa==0 || length(current_param)==0

            elseif len_pa==length(current_param)
                global_phase=subs(global_phase,symvar(global_phase),current_param);   %Change matrix 
            else
                error(['Not the right amount of parameter inputs for global phase.']);
            end
        else %Is a cell
            if mod(length(current_param),2)
                error(['Not the right amount of parameter inputs for global phase. Has to be divisible by 2.'])
            else
                sp=current_param{2:2:end};
                p=current_param{1:2:end};
                global_phase=subs(global_phase,p,sp);   %Change matrix
            end       
        end
        if length(symvar(global_phase))==0
            global_phase=double(global_phase);
        end
        
        if iscell(current_param)==0
            anc_matrix=Embed_Gate(subs(matrix{1},current_param),index{1},N)*global_phase;
        else
            anc_matrix=Embed_Gate(subs(matrix{1},p,sp),index{1},N)*global_phase;
        end
    else
        if iscell(current_param)==0 %Is not a cell -> Use variables as input
            anc_matrix=Embed_Gate(matrix{1},index{1},N);
        else
            if iscell(current_param)==0 
                anc_matrix=Embed_Gate(subs(matrix{1},current_param),index{1},N);
            else
                sp=current_param{2:2:end};
                p=current_param{1:2:end};
                anc_matrix=Embed_Gate(subs(matrix{1},p,sp),index{1},N);
            end
        end
    end
    
    if iscell(current_param)==0
        nums=floor(length(matrix)/fac);
        lens=length(matrix);
        if nums<1
            for i=2:length(matrix)  
                anc_matrix=Embed_Gate(subs(matrix{i},current_param),index{i},N)*anc_matrix;
            end
        else
            timer=tic();
            for i=2:fac
                if isa(matrix{i},'function_handle')==0
                    anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
                else
                    anc_matrix=Embed_Gate(subs(matrix{i},current_param),index{i},N)*anc_matrix;
                end
            end
            t=toc(timer);
            if plotting==1
                fprintf(['    -> Remaining time: ' num2str(lens/fac*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac) '/' num2str(lens) '\n']);
            end
            for j=2:nums
                for i=fac*(j-1)+1:fac*j
                    if isa(matrix{i},'function_handle')==0
                        anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
                    else
                        anc_matrix=Embed_Gate(subs(matrix{i},current_param),index{i},N)*anc_matrix;
                    end
                end
                t=toc(timer);
                if plotting==1
                    fprintf(['    -> Remaining time: ' num2str(lens/(fac*j)*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac*j) '/' num2str(lens) '\n']);
                end
            end
            for i=fac*nums+1:length(matrix)
                if isa(matrix{i},'function_handle')==0
                    anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
                else
                    anc_matrix=Embed_Gate(subs(matrix{i},current_param),index{i},N)*anc_matrix;
                end
            end
            t=toc(timer);
            if plotting==1
                fprintf(['    -> Remaining time: 0s, Elapsed: ' num2str(t) 's, ' num2str(lens) '/' num2str(lens) '\n']);
            end
        end
    else
        nums=floor(length(matrix)/fac);
        lens=length(matrix);
        if nums<1
            for i=2:length(matrix)  
                anc_matrix=Embed_Gate(subs(matrix{i},p,sp),index{i},N)*anc_matrix;
            end
        else
            timer=tic();
            for i=2:fac
                if isa(matrix{i},'function_handle')==0
                    anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
                else
                    anc_matrix=Embed_Gate(subs(matrix{i},p,sp),index{i},N)*anc_matrix;
                end
            end
            t=toc(timer);
            if plotting==1
                fprintf(['    -> Remaining time: ' num2str(lens/fac*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac) '/' num2str(lens) '\n']);
            end
            for j=2:nums
                for i=fac*(j-1)+1:fac*j
                    if isa(matrix{i},'function_handle')==0
                        anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
                    else
                        anc_matrix=Embed_Gate(subs(matrix{i},p,sp),index{i},N)*anc_matrix;
                    end
                end
                t=toc(timer);
                if plotting==1
                    fprintf(['    -> Remaining time: ' num2str(lens/(fac*j)*t-t) 's, Elapsed: ' num2str(t) 's, ' num2str(fac*j) '/' num2str(lens) '\n']);
                end
            end
            for i=fac*nums+1:length(matrix)
                if isa(matrix{i},'function_handle')==0
                    anc_matrix=Embed_Gate(matrix{i},index{i},N)*anc_matrix;
                else
                    anc_matrix=Embed_Gate(subs(matrix{i},p,sp),index{i},N)*anc_matrix;
                end
            end
            t=toc(timer);
            if plotting==1
                fprintf(['    -> Remaining time: 0s, Elapsed: ' num2str(t) 's, ' num2str(lens) '/' num2str(lens) '\n']);
            end
        end
    end
end

ind=1:2^sizes;
anti_ind=2^sizes+1:(2^anc_sizes);
matrix=anc_matrix(ind,ind);

anc_mat1=anc_matrix(anti_ind,ind);
anc_mat2=anc_matrix(ind,anti_ind);
%anc_mat3=anc_matrix(anti_ind,anti_ind); %-> Could be checked aswell -> maybe in a later version
if any(max(anc_mat1(:))>10^-6) || any(max(anc_mat2(:))>10^-6) 
    fprintf(['Warning: The ancilla qubits change after the operation (' num2str(max(anc_mat1(:))) ', ' num2str(max(anc_mat2(:))) ')!\n'])
end
    
end
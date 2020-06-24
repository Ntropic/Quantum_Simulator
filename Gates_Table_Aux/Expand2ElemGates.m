function [ gates_table ] = Expand2ElemGates( i,gates_table,elem_gates )
%Expand to elem_gates -> all gates references have been defined in a
%previous hierarchy and can just be copy pasted out -> this is the
%reason for introducing a hierarchy in the first place (allthough you
%will also see modest speed improvements due to the reordering

gate_ind=gates_table(i).steps.gate_ind; %index of gate in [elem,comp] gate lists
gate_name=gates_table(i).steps.gates;   %name of gate corresponding to gate_ind
indi=gates_table(i).steps.index; %qubit indexes of subgates

gate_len=gates_table(i).single_gates+gates_table(i).double_gates;
gates_table(i).exp_steps.index=cell(gate_len,1); 
gates_table(i).exp_steps.param=cell(gate_len,1); 
gates_table(i).exp_steps.matrix=cell(gate_len,1); 
gates_table(i).exp_steps.gate_ind=cell(gate_len,1); 
gates_table(i).exp_steps.gate_names=cell(gate_len,1); 
gates_table(i).exp_steps.gate_twirling=cell(gate_len,1); 
gates_table(i).exp_steps.depth=cell(gate_len,1); 
gates_table(i).exp_steps.progress=0;
for j=1:length(gate_ind)
    %% Comp Gate
    if gate_ind{j}(1)==0;    % comp_gate in gates_table - copy the whole list of subgates
        %Find elementary gates of this comp_gate -> do same substitution for all of them
        g_i=gate_ind{j}(2);
        current_indexes=gates_table(g_i).exp_steps.index;   %gates_info of subgates elem_gates list
        current_matrices=gates_table(g_i).exp_steps.matrix;
        current_parames=gates_table(g_i).exp_steps.param;  
        current_gate_ind=gates_table(g_i).exp_steps.gate_ind;
        current_gate_name=gates_table(g_i).exp_steps.gate_names;
        current_gate_twirl=gates_table(g_i).exp_steps.gate_twirling;
        current_gate_depth=gates_table(g_i).exp_steps.depth;
        
        current_param=gates_table(i).steps.param{j};
        
        if length(gates_table(g_i).exp_steps.global_phase)>0
            if length(gates_table(i).exp_steps.global_phase)==0
                gates_table(i).exp_steps.global_phase=1;
            end
            %Substitute global phase aswell!!! 
            da_phase=gates_table(g_i).exp_steps.global_phase;
            
            if iscell(current_param)==0 %Is not a cell -> Use variables as input
               	if isa(da_phase,'sym')==1
                    pa=strsplit(char(symvar(da_phase)),{'[',']'});
                    if length(pa)>1
                        pa=pa{2};
                        pa=strsplit(pa,', ');
                    end
                end
                len_pa=length(pa); %Current global_phase parameters
                if len_pa==0 || length(current_param)==0
                    
                elseif len_pa==length(current_param)
                	da_phase=subs(da_phase,symvar(da_phase),current_param);   %Change matrix 
                else
                    error(['Not the right amount of parameter inputs for global phase.']);
                end
            else %Is a cell
                if mod(length(current_param),2)
                    error(['Not the right amount of parameter inputs for global phase. Has to be divisible by 2.'])
                else
                    sp=[current_param{2:2:end}];
                    p=[current_param{1:2:end}];
                    da_phase=subs(da_phase,p,sp);   %Change matrix
                end       
            end
            if length(symvar(da_phase))==0
                da_phase=double(da_phase);
            end
            gates_table(i).exp_steps.global_phase=gates_table(i).exp_steps.global_phase*da_phase;
        end
        
        for k=1:length(current_indexes) %Iterate through all elem_gates of subgate
        	index=indi{j}(current_indexes{k});%New Indexes
            fun_vars=current_parames{k};
            fun_mat=current_matrices{k};
            fun_ind=current_gate_ind{k};
            fun_name=current_gate_name{k};
            fun_twirl=current_gate_twirl{k};
            fun_depth=current_gate_depth{k};
            fun_depth=[j fun_depth(:)'];
            if iscell(current_param)==0 %Is not a cell -> Use variables as input
                if length(current_param)>0
                    len_pa=length(current_parames{k});
                    if len_pa==0 %No substitutions
                        gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,fun_ind,fun_name,fun_twirl,fun_depth);
                    elseif len_pa==length(current_param) %Substitutions without declared variable -> using alphabetic ordering
                        %non cell substitutions are not recommended for
                        %gates with hierarchy>2 because it is not clear
                        %which variables are present in a gate
                        [fun_mat fun_vars]=fun_mat_subs(fun_mat,fun_vars,current_param);   %Change matrix
                        gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,fun_ind,fun_name,fun_twirl,fun_depth);
                    else
                        error(['Not the right amount of parameter inputs for ' gates_table(i).names 's subgate number ' num2str(j) 's elem_gate number ' num2str(k) '.']);
                    end
                else
                    gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,fun_ind,fun_name,fun_twirl,fun_depth);
                end
            else %Is a cell
                if mod(length(current_param),2)
                    error(['Not the right amount of parameter inputs for ' char(gates(i)) '. Has to be divisible by 2.'])
                else
                    [fun_mat fun_vars]=fun_mat_subs(fun_mat,fun_vars,current_param);   %Substitute multiple variables at once
                    gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,fun_ind,fun_name,fun_twirl,fun_depth);
                end           
            end
        end
    %% Elem Gate
    else                    % elem_gate - just add this gate
        index=indi{j};  %Hope this is the correct size... in your gate definitions!
        %Determine substitutions for the elementary gate
        current_param=gates_table(i).steps.param{j};
        ind=gates_table(i).steps.gate_ind{j}(1);
        len_pa=length(elem_gates(ind).fun_vars);
        fun_vars=elem_gates(ind).fun_vars;
        fun_mat=elem_gates(ind).fun_mat;
        fun_twirl=elem_gates(ind).twirling;
        fun_depth=[j];
        if isa(elem_gates(ind).names,'cell')
            name=elem_gates(ind).names{1};
        else
            name=elem_gates(ind).names;
        end
        if iscell(current_param)==0 %Is not a cell -> Use variables as input   
            if len_pa==0 || length(current_param)==0 %No substitutions
                gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,ind,name,fun_twirl,fun_depth);
            elseif len_pa==length(current_param) %Substitutions without declared variable -> using alphabetic ordering
                [fun_mat fun_vars]=fun_mat_subs(fun_mat,fun_vars,current_param);   %Change matrix
                gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,ind,name,fun_twirl,fun_depth);
            else
                error(['Not the right amount of parameter inputs for ' gates_table(i).names ' s subgate ' elem_gates(gates_table(i).exp_steps.index{j}(1)).names '.']);
            end

        else %Is a cell
            if mod(length(current_param),2)
                error(['Not the right amount of parameter inputs for ' char(gates(i)) '. Has to be divisible by 2.'])
            else
                if len_pa~=0 || length(current_param)==0   
                    [fun_mat fun_vars]=fun_mat_subs(fun_mat,fun_vars,current_param);   %Substitute multiple variables at once
                    gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,ind,name,fun_twirl,fun_depth);
                else %No substitutions       
                    gates_table=Exp_Step_Add(i,gates_table,index,fun_vars,fun_mat,ind,name,fun_twirl,fun_depth);
                end
            end           
        end 
    end
end
end


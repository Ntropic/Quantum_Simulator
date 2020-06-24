function [ phi ] = Gate2Matrix_Return_Wave( phi0, elem_gates, other_gates, new_gate, expand_how_deep )
%GATE2MATRIX_Return_Wave creates a matrix from elementary gates 
%(elem_gates), composite and returns the inttermediary pure states of an
%initial phi0 getting transformed by the gate(s)
% Input: 
%   elem_gates, comp_gates  ->  gates tables for elementary and composite 
%                               gates
%	new_gate                ->  gate structure for gate that will be
%                               transformed into matri representation
%   expand_how_deep         ->  How many levels deep should the gates be 
%                               expanded? (-1= Infinity) 
%   re_calc                 ->  (not necessary) 1 recalc anyways
%
% Output:
%   matrix                  ->  The new matrix without ancilla qubits
%   anc_matrix              ->  The new matrix including ancilla qubits
%   gates_expanded          ->  List of gates expanded to depth specified
%                               by (expanded_how deep)
%   gates_indexes_expanded  ->  Indexes for the prev. mentioned
%                               gates_expanded
%   global_phase            ->  product of global phases of the expanded
%                               gates
if nargin==4
    expand_how_deep=1;
end

sizes=new_gate.size;
anc_sizes=new_gate.anc_size;
index=[new_gate.steps.index];
gates=[new_gate.steps.gates];
% Global phase will be multiplied to gate, in order to get a prettier
% matrix representation
if any(strcmp(fields(new_gate),'global_phase'))==1
    if length(new_gate.global_phase)==1
        global_phase=new_gate.global_phase;
    else
        global_phase=1;
    end
else
    global_phase=1;
end
        
%Search the gates in the list
gates_expanded={};
gates_matrices={};
gates_indexes_expanded={};


    for i=1:length(gates)
        %% Is it an elementary gate?
        found=0;
        for j=1:length(elem_gates)
            obj=strcmpi(elem_gates(j).names,gates(i)); %i for ignoring case
            if iscell(obj)
                obj=cell2mat(obj);
            end
            if any(obj)
                gates_expanded={gates_expanded{:},gates(i)};
                gates_indexes_expanded={gates_indexes_expanded{:},[new_gate.steps.index{i}]};
                if length(elem_gates(j).matrix)>0
                    current_param=new_gate.steps.param{i};
                    if length(current_param)~=0
                        if isa(elem_gates(j).matrix,'sym')  %Is symbol or symbolic function
                            %Check if parameters are in string| cell form {'name',factor}
                            if iscell(current_param)==0 %Is not a cell -> Use variables as input
                                elem_gates_curr_matrix=elem_gates(j).matrix;
                                if length(symvar(elem_gates(j).matrix))==0
                                    gates_matrices={gates_matrices{:},elem_gates(j).matrix};
                                elseif length(symvar(elem_gates(j).matrix))==length(current_param)
                                    gates_matrices={gates_matrices{:},subs(elem_gates_curr_matrix,symvar(elem_gates_curr_matrix),current_param)};%elem_gates_curr_matrix(current_param(:))};
                                else
                                    error(['Not the right amount of parameter inputs for ' char(gates(i)) '.']);
                                end
                            else %Is a cell
                                if mod(length(current_param),2)
                                    error(['Not the right amount of parameter inputs for ' char(gates(i)) '. Has to be divisible by 2.'])
                                else
                                    elem_gates_curr_matrix=elem_gates(j).matrix;
                                    k=length(current_param)/2;
                                    symbol_name={current_param{2*k-1}};
                                    form_value={current_param{2*k}};
                                    gates_matrices={gates_matrices{:},subs(elem_gates_curr_matrix,symbol_name,form_value)};
                                end
                            end
                        else %is classical matrix
                            gates_matrices={gates_matrices{:},elem_gates(j).matrix;};
                        end
                    else
                        gates_matrices={gates_matrices{:},elem_gates(j).matrix};
                    end
                else
                    error(['Matrix for elementary gate ' char(elem_gates(j).names) ' not found.'])
                end
            found=1;
            break
            end
        end
        if found==0
            %% Is it a composite gate
            for j=1:length(other_gates)
                obj=strcmpi(other_gates(j).names,gates(i)); %i for ignoring case
                if iscell(obj)
                    obj=cell2mat(obj);
                end
                if any(obj)
                    %Check if matrix exists
                    if isempty(other_gates(j).matrix) || expand_how_deep>0
                        %Recursively use Gate2Matrix
                        [other_gates(j).matrix other_gates(j).anc_matrix gates_expanded_rec gates_indexes_expanded_rec ~]=Gate2Matrix(elem_gates,other_gates,other_gates(j),expand_how_deep-1);
                    end

                    if isa(other_gates(j).matrix,'sym')
                        current_param=new_gate.steps.param{i};
                        current_syms=symvar(other_gates(j).matrix);
                        if iscell(current_param)==0
                            %Check the length of current params
                            if length(symvar(other_gates(j).matrix))==0
                                gates_matrices={gates_matrices{:},other_gates(j).matrix};
                            elseif length(symvar(other_gates(j).matrix))==length(current_param)
                                gates_matrices={gates_matrices{:},subs(other_gates(j).matrix,symvar(other_gates(j).matrix),current_param)};%elem_gates_curr_matrix(current_param(:))};
                            elseif length(current_param)==0
                                gates_matrices={gates_matrices{:},other_gates(j).matrix};
                            else
                                error(['Not the right amount of parameter inputs for ' char(gates(i)) '.']);
                            end
                        else %Is cell
                            if mod(length(current_param),2)
                                error(['Not the right amount of parameter inputs for ' char(gates(i)) '. Has to be divisible by 2.'])
                            else
                                k=1:length(current_param)/2;
                                symbol_name={current_param{2*k-1}};
                                form_value={current_param{2*k}};
                                gates_matrices={gates_matrices{:},subs(other_gates(j).matrix,symbol_name,form_value)};
                            end
                        end
                    else %Is not a symbol
                        gates_matrices={gates_matrices{:},other_gates(j).matrix};
                    end

                    if expand_how_deep>0
                        if iscell(gates_expanded_rec)==0
                            gates_expanded={gates_expanded{:},gates_expanded_rec};
                        else
                            gates_expanded={gates_expanded{:},gates_expanded_rec{:}};
                        end
                        gates_indexes_expanded={gates_indexes_expanded{:},gates_indexes_expanded_rec{:}};
                    end

                    found=1;
                    break
                end
            end
        end
        if found==0
            error(['Gate ' cell2mat(gates(i)) ' not found.']);
        end
    end
    
    phi=zeros(length(phi0),length(gates_matrices)+1);
    phi(:,1)=phi0;
    
    N_anc=sizes+anc_sizes;
    anc_matrix=diag(ones(2^(N_anc),1));

    %Calculate phi's
    curr_i=1;
    for i=length(gates_matrices):-1:1        %From right to left 
        phi(:,curr_i+1)=Embed_Gate( gates_matrices{i} ,index{i},N_anc )*phi(:,curr_i);
        curr_i=curr_i+1;
    end

end
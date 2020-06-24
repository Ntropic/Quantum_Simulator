function [gates_matrices gates_expanded gates_indexes_expanded N N_anc] = Gates_Expand( elem_gates, other_gates, new_gate)
%GATES_EXPAND creates a gate expansion into elementary gates (elem_gates)


sizes=new_gate.size;
anc_sizes=new_gate.anc_size;
index=[new_gate.steps.index];
gates=[new_gate.steps.gates];

%Search the gates in the list
gates_matrices={};
matrix_list={};

gates_expanded={};
gates_indexes_expanded={};

    for i=1:length(gates)
        %Is it an elementary gate?
        found=0;
        for j=1:length(elem_gates)
            obj=strcmpi(elem_gates(j).names,gates(i));
            if iscell(obj)
                obj=cell2mat(obj);
            end
            if any(obj)
                gates_expanded={gates_expanded{:},gates(i)};
                gates_indexes_expanded={gates_indexes_expanded{:},[new_gate.steps.index{i}]};
                if length(elem_gates(j).matrix)>0
                    current_param=new_gate.steps.param{i};
                    if isa(elem_gates(j).matrix,'sym')  %Is symbol or symbolic function
                        %Check if parameters are in string| cell form {'name',factor}
                        if iscell(current_param)==0 %Is not a cell -> Use variables as input
                            if length(symvar(elem_gates(j).matrix))==0
                                gates_matrices={gates_matrices{:},elem_gates(j).matrix};
                            elseif length(symvar(elem_gates(j).matrix))==length(current_param)
                                elem_gates_curr_matrix=elem_gates(j).matrix;
                                gates_matrices={gates_matrices{:},subs(elem_gates_curr_matrix,symvar(elem_gates_curr_matrix),current_param)};
                            else
                                error(['Not the right amount of parameter inputs for ' char(gates(i)) '.']);
                            end
                        else %Is a cell
                            if mod(length(current_param),2)
                                error(['Not the right amount of parameter inputs for ' char(gates(i)) '. Has to be divisible by 2.']);
                            else
                                elem_gates_curr_matrix=elem_gates(j).matrix;
                                for k=1:length(current_param)/2
                                    symbol_name=current_param{2*k-1};
                                    form_value=current_param{2*k};
                                    gates_matrices={gates_matrices{:},subs(elem_gates_curr_matrix,symbol_name,form_value)};
                                end
                            end
                        end
                    else %is classical matrix
                        gates_matrices={gates_matrices{:},elem_gates(j).matrix;};
                    end
                else
                    error(['Matrix for elementary gate ' char(elem_gates(j).names) ' not found.'])
                end

                found=1;
                break
            end
        end
        if found==0
            %Search composite gates
            if iscell(other_gates)==0
                comp_gates=other_gates;
                other_gates=cell(1);
                other_gates{1}=comp_gates;
            end 
            for q=1:length(other_gates)
                comp_gates=other_gates{q};
                for j=1:length(comp_gates)
                    obj=strfind(comp_gates(j).names,gates(i));
                    if iscell(obj)
                        obj=cell2mat(obj);
                    end
                    if any(obj)
                        [ gates_matrices_rec gates_expanded_rec gates_indexes_expanded_rec ]=Gates_Expand(elem_gates,comp_gates,comp_gates(j));
                        
                        gates_expanded={gates_expanded{:},gates_expanded_rec{:}};
                        for o=1:length(gates_indexes_expanded_rec)
                            gates_indexes_expanded={gates_indexes_expanded{:},index{i}(gates_indexes_expanded_rec{o})};
                        end

                        %gates_matrices={gates_matrices{:},gates_matrices_rec{:}};
                        for o=1:length(gates_matrices_rec)
                            if isa(gates_matrices_rec{o},'sym')
                                current_param=new_gate.steps.param{i};
                                if iscell(current_param)==0
                                    %Check the length of current params
                                    if length(symvar(gates_matrices_rec{o}))==0
                                        gates_matrices={gates_matrices{:},gates_matrices_rec{o}};
                                    elseif length(symvar(gates_matrices_rec{o}))==length(current_param)
                                        gates_matrices={gates_matrices{:},subs(gates_matrices_rec{o},symvar(gates_matrices_rec{o}),current_param)};
                                    elseif length(current_param)==0
                                        gates_matrices={gates_matrices{:},gates_matrices_rec{o}};
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
                                        gates_matrices={gates_matrices{:},subs(gates_matrices_rec{o},symbol_name,form_value)};
                                    end
                                end
                            else %Is not a symbol
                                gates_matrices={gates_matrices{:},gates_matrices_rec{o}};
                            end
                        end
                            
                        found=1;
                        break
                    end
                end
            end
        end
        if found==0
            error(['Gate ' cell2mat(gates(i)) ' not found.']);
        end
    end

    N=sizes;
    N_anc=sizes+anc_sizes;
    

end
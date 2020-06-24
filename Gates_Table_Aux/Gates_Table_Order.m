function [ gates_table ] = Gates_Table_Order( gates_table,elem_gates )
%GATES_TABLE_ORDER creates a hierarchy based ordering of a gates table and
%if necessary a hierarchy. (it also updates hierarchies if gates have been
%added)

if any(strcmp(fields(gates_table),'hierarchy'))==0
    [gates_table(:).hierarchy]=deal([]);
end

%First find index representation
for i=1:length(gates_table)
    if any(strcmp(fields(gates_table(i)),'steps'))==1
        %Go through every gate and find index within other_gates,
        for j=1:length(gates_table(i).steps.gates)
            gate_name=gates_table(i).steps.gates{j};
            index=Gate_Index_by_Name(gate_name,elem_gates,gates_table);
            gates_table(i).steps.gate_ind{j}=index;
        end
        if length(gates_table(i).steps.gates)==0
            gates_table(i).steps.gate_ind={};
        end
    end
    gates_table(i).exp_steps.global_phase=gates_table(i).global_phase;
end

%Create gate hierarchy
for i=1:length(gates_table)
    gates_table=Hierarchy_Finder(i,elem_gates,gates_table,1);
end

%Reorder struct by order
T=struct2table(gates_table);
T=sortrows(T,'hierarchy');
gates_table=table2struct(T)';


%Redo the index representation due to new ordering sheme
for i=1:length(gates_table)
    if any(strcmp(fields(gates_table(i)),'steps'))==1
        %Go through every gate and find index within other_gates,
        for j=1:length(gates_table(i).steps.gates)
            gate_name=gates_table(i).steps.gates{j};
            index=Gate_Index_by_Name(gate_name,elem_gates,gates_table);
            gates_table(i).steps.gate_ind{j}=index;
        end
    end
end

end


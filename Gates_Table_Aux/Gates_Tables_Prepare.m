function [ gates_table ] = Gates_Tables_Prepare( gates_table,elem_gates )
%GATES_TABLE_PREPARE prepares gates tables for substitution expansions
%(Gates_Expand_Subs) 

[elem_gates,gates_table]=Gate_Counter(elem_gates,gates_table);

%% Create expanded gates list and substitution parameter lists for each of
%the gates in the expanded gates list
%  Expand gates to elementary gates (using the ordering by hierarchy we 
%  make sure every subgate being called has already been expanded)
i=1;
while i<=length(gates_table)
    gates_table=Expand2ElemGates(i,gates_table,elem_gates);
    i=i+1;
end

end
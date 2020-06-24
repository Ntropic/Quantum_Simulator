function [ elem_gates, gates_table ] = Simple_Errors( single_error,double_error, elem_gates, gates_table )
%SIMPLE_ERRORS creates a simple error model, where all elementary single 
%qubit share the same error rate and elementary  double qubit gates all 
%share the same error rate. combinations of gates get the sum of their 
%constituent error rates. This allows simple cost estimations for substitutions 

gates_table=Gates_Table_Order(gates_table,elem_gates);

for i=1:length(elem_gates(:))
    s=elem_gates(i).size;
    if s==1
        elem_gates(i).cost=single_error;
    else
        elem_gates(i).cost=double_error;
    end
end

for i=1:length(gates_table(:))
    c=0;
    for j=1:length(gates_table(i).steps.gates)
        g=Gate_by_Name(gates_table(i).steps.gates{j},elem_gates,gates_table);
        c=c+g.cost;
    end
    gates_table(i).cost=c;
end

end


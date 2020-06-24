function [ elem_gates, gates_table ] = Gate_Counter( elem_gates, gates_table )
%GATE_COUNTER counts single and double qubit gates within the defined gates

gates_table=Gates_Table_Order(gates_table,elem_gates);

for i=1:length(elem_gates(:))
    s=elem_gates(i).size;
    if s==1
        elem_gates(i).single_gates=1;
        elem_gates(i).double_gates=0;
    else
        elem_gates(i).single_gates=0;
        elem_gates(i).double_gates=1;
    end
end

for i=1:length(gates_table(:))
    c=0;
    d=0;
    for j=1:length(gates_table(i).steps.gates)
        g=Gate_by_Name(gates_table(i).steps.gates{j},elem_gates,gates_table);
        c=c+g.single_gates;
        d=d+g.double_gates;
    end
    gates_table(i).single_gates=c;
    gates_table(i).double_gates=d;
end

end
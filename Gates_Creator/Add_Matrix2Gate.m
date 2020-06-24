function [ gate,subs_gates,comp_gates ] = Add_Matrix2Gate( gate,elem_gates,comp_gates )
%ADD_MATRIX2GATE adds a calculated matrix 2 a gate 

[ matrix anc_matrix subs_gates comp_gates] = Gate2Matrix_Pre_Sub( elem_gates,comp_gates, gate);

if isa(matrix,'sym')
    matrix=simplify(matrix);
end
gate.matrix=matrix;

if isa(anc_matrix,'sym')
    anc_matrix=simplify(anc_matrix);
end
gate.anc_matrix=anc_matrix;

end


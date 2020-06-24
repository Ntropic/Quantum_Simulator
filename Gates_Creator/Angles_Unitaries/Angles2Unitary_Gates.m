function [ U_gates ] = Angles2Unitary_Gates( elem_gates,comp_gates,alpha,beta,delta,theta )
%ANGLES2UNITARY wandelt die Winkel alpha, beta, delta und theta in Unitäre
%2x2 Matrix um. ->Unitary2Angles can calculate the angles from a matrix

    U_gates.size=1;
    U_gates.anc_size=0;
    U_gates.step_num=8;
    U_gates.steps.index={1,1,1,1,1,1,1,1};
    U_gates.steps.gates={'EZ','X','EZ','EY','X','EY','EZ','P01'};
    U_gates.steps.param={[{sym('phi'),1/2*(beta-alpha)}],[],[{sym('phi'),-1/2*(beta+alpha)}],[{sym('phi'),-1/2*theta}],[],[{sym('phi'),1/2*theta}],[{sym('phi'),alpha}],[{sym('phi'),delta}]};

    U_gates.matrix=Gate2Matrix(elem_gates,comp_gates,[],U_gates,2);
end
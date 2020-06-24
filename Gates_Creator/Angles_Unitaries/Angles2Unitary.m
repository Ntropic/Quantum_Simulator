function [ U_gates ] = Angles2Unitary( elem_gates,comp_gates,alpha,beta,delta,theta,name )
%ANGLES2UNITARY wandelt die Winkel alpha, beta, delta und theta in Unitäre
%2x2 Matrix um. 
    if nargin==7
        U_gates.names=name;
    end
    U_gates.size=1;
    U_gates.anc_size=0;
    U_gates.step_num=3;
    U_gates.steps.index={1,1,1};
    U_gates.steps.gates={'EZ','EY','EZ'};
    U_gates.steps.param={[{sym('phi'),beta}],[{sym('phi'),theta}],[{sym('phi'),alpha}]};
    U_gates.global_phase=exp(1i*delta);
    
    U_gates.matrix=Gate2Matrix(elem_gates,comp_gates,U_gates,2);
end
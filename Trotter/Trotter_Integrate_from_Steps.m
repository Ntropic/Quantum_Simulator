function [ A ] = Trotter_Integrate_from_Steps( t,l,step_gate )
%TROTTER_INTEGRATE_FROM_STEPS integrates gates by applying a step_gate l
%times
A=Expansion2Matrix(step_gate,t/l);
A=A^l;
end
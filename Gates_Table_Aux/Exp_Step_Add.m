function [ gates_table ] = Exp_Step_Add( i,gates_table,index,fun_vars,fun_mat,ind,name,twirl,depth )
%Add a step to the exp_step expansion into elementary gates

a=gates_table(i).exp_steps.progress+1;
gates_table(i).exp_steps.progress=a;

gates_table(i).exp_steps.index{a}=index;
if length(fun_vars)==0
    gates_table(i).exp_steps.param{a}=[];
else
    gates_table(i).exp_steps.param{a}=fun_vars;
end
gates_table(i).exp_steps.matrix{a}=fun_mat;

%Also add the gate names and index
gates_table(i).exp_steps.gate_ind{a}=ind;
gates_table(i).exp_steps.gate_names{a}=name;
gates_table(i).exp_steps.gate_twirling{a}=twirl;
gates_table(i).exp_steps.depth{a}=depth;
end


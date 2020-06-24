function [ fun ] = str2fun_mat( c_mat,c_vars )
%STR2FUN_MAT Creates a matlabFunction from a cell of strings c_mat and
%variables in c_vars
c_size=size(c_mat);

c_mat=reshape(c_mat,[1,c_size(1)*c_size(2)]);
c_mat=strjoin(c_mat,',');

c=[c_vars,'[',c_mat,'],[',num2str(c_size(1)),',',num2str(c_size(2)),'])'];

fun=str2func(c);


end
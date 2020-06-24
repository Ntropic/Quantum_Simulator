function [ c_mat,c_vars ] = fun_mat2str( fun )
%FUN_MAT2STR Creates a cell of strings of a matlabFunction containing a
%matrix

c=char(fun);
c=strsplit(c,{'[',']'});
c_vars=c{1};
c_mat=c{2};
c_size=str2num(['[' c{4} ']']);

c_mat=strsplit(c_mat,',');
c_mat=reshape(c_mat,c_size);

end
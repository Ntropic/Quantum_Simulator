function FockSparsePrint( S,extras,extras2,extras3 )
%FOCKSPARSEPRINT prints the matrix as a fock matrix with default 
%substitutions for unitary rotations
if nargin==1
    [~,~,M_mode]=Sparse2Square(S,'fock_unit_def_subs');
    fprintf(M_mode)
elseif nargin==2
    [~,~,M_mode]=Sparse2Square(S,['fock_unit_def_subs_' extras]);
    fprintf(M_mode)
elseif nargin==3
    %extras2(1) max number of photons per mode n
    %extras2(2) number of modes N
    [~,~,M_mode]=Sparse2Square(S,{['unit_def_subs_fock_' extras],extras2(1),extras2(2)});
    fprintf(M_mode)
elseif nargin==4
    %extras2(1) max number of photons per mode n
    %extras2(2) number of modes N
    [~,~,M_mode]=Sparse2Square(S,{['unit_def_subs_fock_' extras],extras2(1),extras2(2)},'index_list',extras3);
    fprintf(M_mode)
end
end


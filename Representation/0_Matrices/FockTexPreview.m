function FockTexPreview( S,name,extras,extras2,extras3 )
%FOCKPRINT prints the matrix as a fock matrix with default substitutions for
%unitary rotations

if nargin==1
    [f]=Sparse2Tex(S,'test','tikz_fock_def_subs');
    lualatex2pdf2preview(f);
elseif nargin==2
    [f]=Sparse2Tex(S,name,'tikz_fock_def_subs');
    lualatex2pdf2preview(f);
elseif nargin==3
    [f]=Sparse2Tex(S,name,['tikz_fock_def_subs_' extras]);
    lualatex2pdf2preview(f);
elseif nargin==4
    %extras2(1) max number of photons per mode n
    %extras2(2) number of modes N
    [f]=Sparse2Tex(S,name,{['tikz_def_subs_fock_' extras],extras2(1),extras2(2)});
    lualatex2pdf2preview(f);
elseif nargin==5
    %extras2(1) max number of photons per mode n
    %extras2(2) number of modes N
    [f]=Sparse2Tex(S,name,{['tikz_def_subs_fock_' extras],extras2(1),extras2(2)},'index_list',extras3);
    lualatex2pdf2preview(f);
end
end


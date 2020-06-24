function [ file_names ] = Gates2Circ2Preview( gate,elem_gates,comp_gates,file_name,mode,depth,sep_dist )
%Creates a qcircuit LaTeX file and compiles it, creating a picture and
%previewing it
%mode specifies the kind of gate that will be constructed
%       - 'gate' creates the gate
%       - 'expansion' creates the expansion of the gate
%       - 'no_sep' the structure does not get separated into multiple plots
%       - 'no_head' delete the header of the tex file
%       - gate=expansion' creates both next to each other (not yet implemented)
%               Can be done via ghost qubits for missing ancilla bits
%depth is how deep the hierarchies get unraveled
%sep_dist is how many gate columns are shown per file
if nargin<4
    error('Not enough input arguments.')
end
if nargin==4
    mode='gate=expansion';
    depth=1;
    sep_dist=20;
elseif nargin==5
    depth=1;
    sep_dist=20;
elseif nargin==6
    sep_dist=20;
end
    
file_names=Gates2Tex(gate,elem_gates,comp_gates,file_name,mode,depth,sep_dist);

%Plot the results
for i=1:length(file_names)
    tex2pdf2preview(file_names{i},1,400);
end
end


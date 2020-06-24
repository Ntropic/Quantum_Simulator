function [ conc_gates ] = Comp_Gate_Merger( comp_gates, comp_gates2,varargin )
%COMP_GATES_MERGER merges/concatenates different compositional gates tables 
%using only the fields shared between the first entered gates table and the second. 

%Find field names of table comp_gates 
fn=fieldnames(comp_gates)';
if length(comp_gates2)>0
    fn2=fieldnames(comp_gates2)';
else
    fn2=[];
end
[a,b,c]=intersect(fn,fn2);
bs=sort(b);

len=length(comp_gates);
    
comp_gates=comp_gates(:);
comp_gates2=comp_gates2(:);
    
conc_gates=struct();
for i=fn(bs)
    for l=1:length(comp_gates)
       conc_gates(l).(i{1})=comp_gates(l).(i{1});
    end
    for l=1:length(comp_gates2)
        conc_gates(l+len).(i{1})=comp_gates2(l).(i{1});
    end
end

%Also add field names of non intersecting fields
[d,e]=setdiff(fn,fn2);
for i=fn(e)
    for l=1:length(comp_gates)
       conc_gates(l).(i{1})=comp_gates(l).(i{1});
    end
end

[d,e]=setdiff(fn2,fn);
for i=fn2(e)
    for l=1:length(comp_gates2)
       conc_gates(l+len).(i{1})=comp_gates2(l).(i{1});
    end
end


%Recursion
if nargin>3
    conc_gates=Comp_Gate_Merger( conc_gates, varargin{1},varargin{2:end});
elseif nargin==3
    conc_gates=Comp_Gate_Merger( conc_gates, varargin{1});
end

end


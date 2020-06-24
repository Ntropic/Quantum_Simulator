function [ s_q_pred ] = Difference_Predictor ( s_q,dim )
%DIFFERENCE_PREDICTOR predicts the gate number increase for gates s_q along
%the Trotter dimension (dim) via difference and averaging

if isa(s_q,'cell')
    sizes=size(s_q{1});
    sizes2=sizes;
    sizes2(dim)=[];
    for i=1:length(s_q)
        d=diff(s_q{i},1,dim);
        s=sum(d,dim)/size(d,dim);
        s_q_pred{i}=reshape(s,sizes2);
    end
else
        d=diff(s_q,1,dim);
        s=sum(d,dim)/size(d,dim);
        s_q_pred=reshape(s,sizes2);
end
end


function [ seq,d_seq ] = Higher_Trotter_Time_Steps( steps )
%HIGHER_TROTTER_TIME_STEPS creates a time stepping map
% 0     -> 3 steps
% 1     -> 5 steps
% n     -> 2n+3 steps
% n+0.5 -> 2n+4 steps
n=1:length(steps);
m=2*n+1;
%Coefficients for different types of recursive approach
s=1./(2-2.^(1./m));
s2=1./(4-4.^(1./m));

seq=[0 1];
for i=1:length(steps)
    if steps(i)==0 %3 substeps
        seq_new=[seq*s(i)];
        seq_new=[seq_new seq_new(end)+(1-2*s(i))*seq(2:end)];
        seq_new=[seq_new seq_new(end)+seq(2:end)*s(i)];    
        
        seq=seq_new;
    elseif steps(i)==1 %5 substeps
        seq_new=[seq*s2(i)];
        seq_new=[seq_new seq_new(end)+seq(2:end)*s2(i)];
        seq_new=[seq_new seq_new(end)+(1-4*s2(i))*seq(2:end)];
        seq_new=[seq_new seq_new(end)+seq(2:end)*s2(i)];
        seq_new=[seq_new seq_new(end)+seq(2:end)*s2(i)];
        
        seq=seq_new;
    elseif mod(steps(i),1)==0 % more substeps -> Generalization of septupling and beyond
        k=2*(steps(i)+1);
        sn=1./(k-k.^(1./m));
        seq_new=[seq*sn(i)];
        for j=1:steps(i)
            seq_new=[seq_new seq_new(end)+seq(2:end)*sn(i)];
        end
        seq_new=[seq_new seq_new(end)+(1-k*sn(i))*seq(2:end)];
        for j=1:steps(i)+1
            seq_new=[seq_new seq_new(end)+seq(2:end)*sn(i)];
        end
        
        seq=seq_new;
    else % more substeps -> Generalization of quartupling and beyond
        steps(i)=floor(steps(i));
        k=2*steps(i)+3;
        sn=1./(k-k.^(1./m));
        seq_new=[seq*sn(i)];
        for j=1:steps(i)
            seq_new=[seq_new seq_new(end)+seq(2:end)*sn(i)];
        end
        seq_new=[seq_new seq_new(end)+(1-k*sn(i))*seq(2:end)];
        for j=1:steps(i)+1
            seq_new=[seq_new seq_new(end)+seq(2:end)*sn(i)];
        end
        
        seq=seq_new;
    end
end
d_seq=diff(seq);
end


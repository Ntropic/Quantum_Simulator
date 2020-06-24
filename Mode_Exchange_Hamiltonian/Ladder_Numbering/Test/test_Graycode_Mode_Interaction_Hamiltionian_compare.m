%test_Graycode_Mode_Interaction_Hamiltionian_compare.m

clear all;
close all;
clc;

do_print=0;
do_tex=0;
do_anim=1;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Number of qubits per mode
Ni=1:3;
ni=2*Ni;
sizle=2.^ni;


for i=1:length(Ni)
    N=Ni(i);
    mode={'fock_gray_gif',2^N-1,2};
    %Create/Load Controlled Unitaries for different sizes
    H=Gray_Exchange_Hamiltonian(N);

    
    h=subplot(1,3,i);
    g=pcolor([H zeros(length(H),1);zeros(1,length(H)+1)]);
    if max(size(H))>64
        set(g, 'EdgeColor', 'none')
    end
    set(h,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis equal;
    axis tight;
    set(h,'XTickLabelRotation',90)
    
    bra_str=bra_fock_str(1:size(H),mode);
    ket_str=ket_fock_str(1:size(H),mode);
    
    max_ticks=16;
    tick_num=min([length(H),max_ticks]);
    
    if tick_num~=length(H) 
        div=length(H)/tick_num; 
        indexes=1:div:length(H); 
    else
        indexes=1:length(H);
    end
    
    if tick_num==length(H) %Take all ticks
        h.XTick=indexes+0.5;
        h.YTick=indexes+0.5;
        h.XTickLabel=ket_str(indexes);
        h.YTickLabel=bra_str(indexes);
    else 
        h.XTick=indexes+0.5;
        h.YTick=indexes+0.5;
        h.XTickLabel=ket_str(indexes);
        h.YTickLabel=bra_str(indexes);
    end
        
end
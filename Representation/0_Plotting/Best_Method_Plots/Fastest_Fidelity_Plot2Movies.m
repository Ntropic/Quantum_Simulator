function  Fastest_Fidelity_Plot2Movies( mode,fid_list,data )
%FASTEST_FIDELITY_PLOT2SETS Plots the fastest Fidelity comparing different
%approaches for 2 datasets
%data={s_q_pred,s_q_pred2,F_pred,F_pred2,F,F2,s_q,s_q2,costs,names1,names2,N_m,N_m2,n_m,n_m2};

  
    h=figure('Name','Fastest Fidelity Plot','Position',[100 100 800 600]);
    ax=subplot(1,1,1);  
    set(gcf,'color',0.98*[1 1 1])

vidfile = VideoWriter('fastest_fidelity.avi','Uncompressed AVI');
vidfile.FrameRate = 5;
open(vidfile);
for i=1:length(fid_list)
    fid=fid_list(i);
    %Extract data
    s_q_pred=data{1};
    s_q_pred2=data{2};
    F_pred=data{3};
    F_pred2=data{4};
    F=data{5};
    F2=data{6};
    s_q=data{7};
    s_q2=data{8};
    costs=data{9};
    names1=data{10};
    names2=data{11};
    N_m=data{12};
    N_m2=data{13};
    n_m=data{14};
    n_m2=data{15};
    cost_maps=data{16};
    cost_maps2=data{17};
    extrap_maps=data{18};
    extrap_maps2=data{19};

    ax=subplot(1,1,1);    
    [cost_maps extrap_maps]=Fastest_Fidelity(fid,F,s_q,costs,F_pred,s_q_pred);
    [cost_maps2 extrap_maps2]=Fastest_Fidelity(fid,F2,s_q2,costs,F_pred2,s_q_pred2);

    data{16}=cost_maps;
    data{17}=cost_maps2;
    data{18}=extrap_maps;
    data{19}=extrap_maps2;

    h=Draw_Fastest_Fidelity_Plot2Sets([h ax], mode,fid,data );
    
    A=getframe(gcf); 
    writeVideo(vidfile,A);
end
close(vidfile)
end
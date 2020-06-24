function h=Compare_Fidelity_Distribution_Circ_XY_All_Interaction()
    %Compare_Fidelity_Distribution_Circ_XY_All_Interaction.m
    clear all;
    close all;
    clc;

    p=pwd;
    if any(strfind(p,'\'));
        elem=strsplit(p,'\');
    else
        elem=strsplit(p,'/');
    end
    shortened=fullfile(elem{1:end-3});
    addpath(genpath(shortened));
    %Load gates
    load('../../../Gates_Table/elem_gates.mat','-mat')
    load('../../../Gates_Table/comp_gates.mat','-mat')
    load('../XY_quantum_gates.mat','-mat')

    %Trotter step number
    steps=20;

    %Iteration time
    t=pi/4;

    n_m=3:10;
    
    
    %Higher Trotter recursion order
    rec_order=1;

    bins=200;
    how_many=10^5;

    for i=1:length(n_m)
        n=n_m(i);
        fprintf([num2str(i) ' Photons:\n'])
        N=n-1;%Photon number

        Hij=XY_Exchange_Hamiltonian(n);
        U_exact=expm(-1i*t*Hij);

        %Subspace for analysis of algorithm
        index=2.^(0:n-1)+1;

        Ni=N:-1:1;
        in=1:N;
        Jx=sqrt(in).*sqrt(Ni);

        %Gate A
        sizle=floor(n/2)*2;
        anc_size=0;
        nameA=['A_n_' num2str(n)];
        circ_string=['A_{XY}^{n=' num2str(n) '}'];
        gates_indexes={};
        param_steps={};
        gates_names={};
        for j=1:2:2*floor(n/2)-1
            gates_indexes={gates_indexes{:},[j,j+1]};
            param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(j)}]};
        end
        [gates_names{1:length(gates_indexes)}]=deal('XX_YY');
        A=Create_Comp_Gate(nameA,sizle,anc_size,gates_names,gates_indexes,param_steps);
        A=Generate_Gate_Circuit(A,circ_string,1:sizle,[],[]);

        %Gate B
        sizle=floor((n-1)/2)*2;
        anc_size=0;
        nameB=['B_n_' num2str(n)];
        circ_string=['B_{XY}^{n=' num2str(n) '}'];
        gates_indexes={};
        param_steps={};
        gates_names={};
        for j=1:2:2*floor((n-1)/2)-1
            gates_indexes={gates_indexes{:},[j,j+1]};
            param_steps={param_steps{:},[{sym('phi'),sym('phi')*Jx(j+1)}]};
        end
        [gates_names{1:length(gates_indexes)}]=deal('XX_YY');
        B=Create_Comp_Gate(nameB,sizle,anc_size,gates_names,gates_indexes,param_steps);
        B=Generate_Gate_Circuit(B,circ_string,1:sizle,[],[]);


        for l=1:steps
            fprintf([' - ' num2str(l) ' Trotter steps\n'])
            timer=tic();

            %Create Trotter steps
            name=['XY_n_' num2str(n)];
            circ_string=['XY_{\text{approx}}^{n=' num2str(n) '}'];
            circ_string_step=['XY_{\text{step}}^{n=' num2str(n) '}'];
            gates_names={nameA,nameB};
            gates_indexes={[1:A.size],1+[1:B.size]};
            other_gates1=Comp_Gate_Merger(comp_gates,xy_gates,A,B);
            %This command creates a Trotter gate for all the steps -> Great for
            %creating a gate circuit
             pattern=rec_order;
            [seq,d_seq,gate,step_gate]=Higher_Trotter(t,l,pattern,'',circ_string,circ_string_step,n,anc_size,gates_names,gates_indexes);
        %     if l==1
        %         Gates2Circ2Preview(gate,elem_gates,other_gates1,'xy_gate','gate=expansion',2,30);
        %     end

            %Do Trotter steps
            other_gates=Comp_Gate_Merger(other_gates1,gate,step_gate);
            other_gates=Gates_Tables_Prepare(other_gates,elem_gates);

            ind=Gate_Index_by_Name(gate.names,elem_gates,other_gates);
            U_approx=Expansion2Matrix(other_gates(ind(2)));

            [F_min(l) F_av(l)]=Min_Av_Fidelity(U_approx,U_exact,index);
            U_ap{l}=U_approx;
            %average_fidelity
            time_l(l)=toc(timer);
        end 

        phi_noon=NOON_Wave_XY(N);
        for l=1:steps
            F_noon(l)=Fidelity(U_ap{l},U_exact,index,phi_noon(index));
        end

        bins_int=[min([1-F_av 1-F_noon]) max([1-F_min 1-F_noon])];
        bins_int=log10([bins_int(2),bins_int(1)]);
        xs=10.^(linspace(bins_int(1),bins_int(2),bins));
        x{i}=xs;

        phi=RandWave(n,how_many);
        for l=1:steps
            F_bins(l,:)=Fidelity(U_ap{l},U_exact,index,phi);

            F_hist(l,:)=hist(1-F_bins(l,:),xs(end:-1:1));
        end
        F_hist=F_hist(:,end:-1:1);


        %Calculate Order of Convergence
        f1_bins{i}=1-F_bins;
        for j=1:how_many
            order2(:,j)=Order_Of_Convergence(f1_bins{i}(:,j));
        end
        order{i}=order2;
        xo_min{i}=min(order{i}(:));
        xo_max{i}=max(order{i}(:));
        bins_o{i}=[xo_min{i} xo_max{i}];
        xo{i}=linspace(bins_o{i}(1),bins_o{i}(2),bins);
        for l=1:steps-1
            or_hi(l,:)=hist(order{i}(l,:),xo{i});
        end
        order_hist{i}=or_hi;

        f1_av{i}=1-F_av;
        f1_min{i}=1-F_min;
        for k=1:steps
            fii(:,k)=F_hist(k,:)'/max(F_hist(k,:));
        end
        f1_bins{i}=fii;
        order_hist{i}=order_hist{i}'/max(order_hist{i}(:));
        f1_noon{i}=1-F_noon;

        Order_f1_av{i}=Order_Of_Convergence(f1_av{i});
        Order_f1_min{i}=Order_Of_Convergence(f1_min{i});
        Order_f1_noon{i}=Order_Of_Convergence(f1_noon{i});
    end

    data={steps,x,xo,Order_f1_av,Order_f1_min,Order_f1_noon,order_hist,f1_bins,f1_av,f1_min,f1_noon,f1_bins,n_m};
    save(['r0_comp_fid_dist_N_' num2str(max(n_m)) '.mat'],'data');
    
    %% Plot the figures
    i=1;
    %% 1-Fidelity Plot
    h=figure('Position',[100 100 1200 800]);
    set(h,'color',0.98*[1 1 1]);

    ax1=subplot(1,2,1);
    lin_x=linspace(min(x{i}),max(x{i}),length(x{i}));
    dx=lin_x(2)-lin_x(1);
    lin_x=lin_x+dx/2-(x{i}(end-1)-x{i}(end));
    imagesc(ax1,1:steps,lin_x,f1_bins{i}(end:-1:1,:),'HandleVisibility','off');
    je=gray(512);
    je=je(255:512,:);
    je=je(end:-1:1,:);
    colormap(je)
    caxis([0 max(max(f1_bins{i}(1:end-1,:)))])
    set(ax1,'YDir','normal')
    set(ax1, 'YScale', 'log')
    hold on;
    plot3(ax1,1:steps,f1_av{i},ones(size(f1_av{i})),'k')
    plot3(ax1,1:steps,f1_min{i},ones(size(f1_min{i})),'--k')
    plot3(ax1,1:steps,f1_noon{i},ones(size(f1_noon{i})),'--r')
    xlabel(ax1,'Trotter steps')
    ylabel(ax1,'1-Fidelity')
    legend(ax1,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax1,[1 steps min(x{i}) max(x{i}) 0 1]);
    view(ax1,0,90)
    hold off;
    drawnow;

    %% Order of Convergence of 1-Fidelity Plot
    ax2=subplot(1,2,2);
    lin_xo=linspace(min(xo{i}),max(xo{i}),length(xo{i}));
    dxo=lin_xo(2)-lin_xo(1);
    lin_xo=lin_xo-dxo/2;
    imagesc(ax2,1:(steps-1),lin_xo(end:-1:1),order_hist{i}(end:-1:1,:),'HandleVisibility','off');

    colormap(je)
    caxis([0 max(max(order_hist{i}))])
    set(ax2,'YDir','normal')
    hold on;
    plot3(ax2,1:(steps-1),Order_f1_av{i},ones(size(Order_f1_av{i})),'k')
    plot3(ax2,1:(steps-1),Order_f1_min{i},ones(size(Order_f1_min{i})),'--k')
    plot3(ax2,1:(steps-1),Order_f1_noon{i},ones(size(Order_f1_noon{i})),'--r')
    xlabel(ax2,'Trotter steps')
    ylabel(ax2,'Order of Convergence')
    legend(ax2,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax2,[1 steps-1 min(xo{i}) max(xo{i}) 0 1]);
    view(ax2,0,90)
    hold off;
    drawnow;
    
    %% UI Elements

    photon_slider=uicontrol('Parent',h,'Style','slider','Tag','ph_slide','Position',[1125 150 20 310],'Min',1,'Max',length(n_m),'Value',1,'SliderStep',[1/length(n_m) 1/length(n_m)]);
    addlistener(photon_slider,'ContinuousValueChange',@(src,evt)Redraw(src,evt,[ax1 ax2],h,data));
    photon_text=uicontrol('Parent',h,'Style','text','Tag','photon_text','Position',[1085 90 100 50],'String',['N=' num2str(n_m(1)-1)],'BackgroundColor',get(h,'Color'),'HorizontalAlignment','center');
    
    
end
function Redraw(a,b,ax,h,data)
    ax1=ax(1);
    ax2=ax(2);
        
    %data={steps,x,xo,Order_f1_av,Order_f1_min,Order_f1_noon,order_hist,f1_bins,f1_av,f1_min,f1_noon,f1_bins};
    steps=data{1};
    x=data{2};
    xo=data{3};
    Order_f1_av=data{4};
    Order_f1_min=data{5};
    Order_f1_noon=data{6};
    order_hist=data{7};
    f1_bins=data{8};
    f1_av=data{9};
    f1_min=data{10};
    f1_noon=data{11};
    f1_bins=data{12};
    n_m=data{13};

    %Find uicontrol elements
    photon_slider=findobj(h,'Tag','ph_slide');
    photon_text=findobj(h,'Tag','photon_text');

    slide_num=round(get(photon_slider,'Value'));  
    set(photon_text,'String',['N=' num2str(n_m(slide_num)-1)]);
    
    i=slide_num;
    %% 1-Fidelity Plot
    lin_x=linspace(min(x{i}),max(x{i}),length(x{i}));
    dx=lin_x(2)-lin_x(1);
    lin_x=lin_x+dx/2-(x{i}(end-1)-x{i}(end));
    hold(ax1,'off');
    caxis([0 max(max(f1_bins{i}(1:end-1,:)))])
    imagesc(ax1,1:steps,lin_x,f1_bins{i}(end:-1:1,:),'HandleVisibility','off');
    je=gray(512);
    je=je(255:512,:);
    je=je(end:-1:1,:);
    colormap(je)

    set(ax1,'YDir','normal')
    set(ax1, 'YScale', 'log')
    hold(ax1,'on');
    
    plot3(ax1,1:steps,f1_av{i},ones(size(f1_av{i})),'k')
    plot3(ax1,1:steps,f1_min{i},ones(size(f1_min{i})),'--k')
    plot3(ax1,1:steps,f1_noon{i},ones(size(f1_noon{i})),'--r')
    xlabel(ax1,'Trotter steps')
    ylabel(ax1,'1-Fidelity')
    %legend(ax1,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax1,[1 steps min(x{i}) max(x{i}) 0 1]);
    view(ax1,0,90)
    hold(ax1,'off');

    %% Order of Convergence of 1-Fidelity Plot
    lin_xo=linspace(min(xo{i}),max(xo{i}),length(xo{i}));
    dxo=lin_xo(2)-lin_xo(1);
    lin_xo=lin_xo-dxo/2;
    hold(ax2,'off');
    caxis([0 max(max(order_hist{i}))])
    imagesc(ax2,1:(steps-1),lin_xo(end:-1:1),order_hist{i}(end:-1:1,:),'HandleVisibility','off');
    colormap(je)

    set(ax2,'YDir','normal')
    hold(ax2,'on');
    Order_f1_av{i}
    plot3(ax2,1:(steps-1),Order_f1_av{i},ones(size(Order_f1_av{i})),'k')
    plot3(ax2,1:(steps-1),Order_f1_min{i},ones(size(Order_f1_min{i})),'--k')
    plot3(ax2,1:(steps-1),Order_f1_noon{i},ones(size(Order_f1_noon{i})),'--r')
    xlabel(ax2,'Trotter steps')
    ylabel(ax2,'Order of Convergence')
    legend(ax2,{'Average Fidelity','Minimum Fidelity','NOON Fidelity'},'Location','northeast')
    axis(ax2,[1 steps-1 min(xo{i}) max(xo{i}) 0 1]);
    view(ax2,0,90)
    hold(ax2,'off');

    %matlab2tikz(['Fidelity_Density/r0_fidelity_density_N_' num2str(N) '.tex'],'standalone',true)
end
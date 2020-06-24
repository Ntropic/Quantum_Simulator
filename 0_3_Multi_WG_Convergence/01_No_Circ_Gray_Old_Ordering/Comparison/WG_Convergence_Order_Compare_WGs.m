%Convergence_Order_Compare_Trotter.m
clear all;
close all;
clc;

mode=2;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Photon number
n_m=1:7;
n_m2=1:3;

%Trotter Iteration number
steps=10;

%Waveguide number
N_m=2:4;
N_m2=4:6; %Loading only depends on maximum

%Iteration time
t=pi/4;


name=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_N_' num2str(max(N_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
name2=['trot_compare_n_' num2str(min(n_m2)) '_to_' num2str(max(n_m2)) '_N_' num2str(max(N_m2)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];

a=load([name '.mat']);
a2=load([name2 '.mat']);

F_av_ind=a.F_av_ind;
F_av_ind2=a2.F_av_ind;

%% Plot Convergence 
% For Waveguides

h_c_wg=figure('Position',[100,100,1400,600]);
c=[];
c2=[];
for l=1:steps
    for j=2:size(F_av_ind,2)
        %(#waveguides , #photons , #Trotter steps)
        c(:,j-1,l)=Order_Of_Convergence(1-F_av_ind(1:end,j,l),N_m);
    end
    for j=2:size(F_av_ind2,2)
        c2(:,j-1,l)=Order_Of_Convergence(1-F_av_ind2(1:end,j,l),2:6);
    end
end

if mode==1
    for i=1:steps;
        ci=reshape(c(:,:,i),size(c,1),size(c,2),1);
        c2i=reshape(c2(2:4,:,i),3,size(c2,2),1);

        surf(2:7,[3 4],ci);
        hold on;
        surf(2:3,N_m2,c2i);

        xlabel('Number of photons')
        zlabel('Convergence Order (# of Waveguides)')
        ylabel('Waveguides')

        ax=get(h_c_wg,'children');
        ax.YTick=[3:6];
        ax.XTick=[2:7];

        view(-45,60)
        axis tight;
        drawnow;
        pause(0.2)
    end
elseif mode==2       
    for i=1:size(c,1)-1
        for j=1:size(c,2)
            surf(j+1+(0:steps-1)*1/9,[i i+1]+2,reshape(c(i:i+1,j,:),2,steps));
            hold on;
        end
    end
    for i=2:size(c2,1)-1
        for j=1:size(c2,2)
            surf(j+1+(0:steps-1)*1/9,[i i+1]+2,reshape(c2(i:i+1,j,:),2,steps));
            hold on;
        end
    end

    xlabel('Number of photons')
    zlabel({'Convergence Order','\fontsize{9}(# of Waveguides)'})
    ylabel('Waveguides')

    ax=get(h_c_wg,'children');
    ax.YTick=[3:6];
    ax.XTick=[1:7];

    z=c(2,4,1)-c(2,4,10);
    
    quiver3(5,4.05,0,1,0,-z,'k')
    
    te=text(5.5,4.15,0,'Trotter steps','HorizontalAlignment','center','FontSize',10);
    set(te,'Rotation',19)
    
    view(-45,60)
    axis tight;
    drawnow;
elseif mode==3      
    ci=reshape(c(:,:,1),size(c,1),size(c,2),1);
    c2i=reshape(c2(2:4,:,1),3,size(c2,2),1);
    cf=reshape(c(:,:,10),size(c,1),size(c,2),1);
    c2f=reshape(c2(2:4,:,10),3,size(c2,2),1);

    ma=mesh(2:7,[3 4],ci,'FaceAlpha',0);
    hold on;
    qa=quiver3(2:7,3:4,ci,zeros(2,6),zeros(2,6),cf-ci,'Color','k');
    
    mb=mesh(2:3,N_m2,c2i,'FaceAlpha',0);
    qb=quiver3(2:3,5:6,c2i(2:3,:),zeros(2,2),zeros(2,2),c2f(2:3,:)-c2i(2:3,:),'Color','k');
    
    xlabel('Number of photons')
    zlabel('Convergence Order (# of Waveguides)')
    ylabel('Waveguides')

    ax=get(h_c_wg,'children');
    ax.YTick=[3:6];
    ax.XTick=[2:7];

    view(-45,60)
    axis tight;
        
    set(qa,'AutoScale','off');
    set(qb,'AutoScale','off');
    set(qa,'MaxHeadSize',0.05);
    set(qb,'MaxHeadSize',0.05);
    drawnow;
elseif mode==4     
    ci=reshape(c(:,:,1),size(c,1),size(c,2),1);
    c2i=reshape(c2(2:4,:,1),3,size(c2,2),1);
    cf=reshape(c(:,:,10),size(c,1),size(c,2),1);
    c2f=reshape(c2(2:4,:,10),3,size(c2,2),1);

    ma=surf(2:7,3:4,ci);
    hold on;
    qa=quiver3(2:7,3:4,ci,zeros(2,6),zeros(2,6),cf-ci,'Color','k');
    maf=surf(2:7,3:4,cf);
    
    mb=surf(2:3,N_m2,c2i);
    qb=quiver3(2:3,5:6,c2i(2:3,:),zeros(2,2),zeros(2,2),c2f(2:3,:)-c2i(2:3,:),'Color','k');
    mbf=surf(2:3,N_m2,c2f);
    
    xlabel('Number of photons')
    zlabel('Convergence Order (# of Waveguides)')
    ylabel('Waveguides')

    ax=get(h_c_wg,'children');
    ax.YTick=[3:6];
    ax.XTick=[2:7];

    view(-45,60)
    axis tight;
        
    set(qa,'AutoScale','off');
    set(qb,'AutoScale','off');
    set(qa,'MaxHeadSize',0.05);
    set(qb,'MaxHeadSize',0.05);
    drawnow;
end
%matlab2tikz('convergence_order_trotter

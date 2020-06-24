%Convergence_Order_Compare_Sym.m
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


name=['sym_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_N_' num2str(max(N_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];
name2=['sym_compare_n_' num2str(min(n_m2)) '_to_' num2str(max(n_m2)) '_N_' num2str(max(N_m2)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];

a=load([name '.mat']);
a2=load([name2 '.mat']);

F_av_ind=a.F_av_ind;
F_av_ind2=a2.F_av_ind;

%% Plot 1-Fidelity
% h=figure('Position',[100,100,1400,600]);
% m=size(F_av_ind,1);
% for i=1:m
%     ax{i}=subplot(1,m,i);
%     err=reshape(1-F_av_ind(i,:,:),size(F_av_ind,2),size(F_av_ind,3));
%     err2=reshape(1-F_av_ind2(i,:,:),size(F_av_ind2,2),size(F_av_ind2,3));
%     surf(err');
%     hold on;
%     surf(err2');
%     if i==1
%         xlabel('Number of photons')
%         zlabel('1-Fidelity')
%     end
%     if i==m
%         ylabel('Trotter steps')
%     end
%     title([num2str(N_m(i)) ' Waveguides']);
%     view(-135,30)
%     axis tight;
%     drawnow;
% end

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

quiver3(5.25,4.15 ,0.1,1,0,0.025,'k')

te=text(5.8,4.275,0,'Trotter steps','HorizontalAlignment','center','FontSize',10);
set(te,'Rotation',19)

view(-45,60)
axis tight;
set(h_c_wg,'color',0.98*[1 1 1])
drawnow();


A=getframe(h_c_wg);
imwrite(A.cdata,'wg_convergence_sym.png');

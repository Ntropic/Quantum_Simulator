%Convergence_Order_Compare_Trotter.m
clear all;
%close all;
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

%Trotter Iteration number
steps=10;

%Waveguide number
N_m=2:4;

%Iteration time
t=pi/4;


name=['trot_compare_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_N_' num2str(max(N_m)) '_s_' num2str(steps) '_t_' num2str(t) '_pauli_gray'];

a=load([name '.mat']);

F_av_ind=a.F_av_ind;

%% Plot Convergence 
% For Photon numbers
h_c_p=figure('Position',[100,100,1400,600]);
m=1;%size(F_av_ind,1);
for i=1%1:m
    ax_h_t{i}=subplot(1,m,i);
    c=[];
    for j=1:size(F_av_ind,3)
        c(j,:)=Order_Of_Convergence(1-F_av_ind(i,2:end,j),2:max(n_m));
    end
    surf(3:max(n_m),1:steps,c);
    if i==1
        xlabel('Number of photons')
        zlabel('Convergence Order (# of Photons)')
    end
    if i==m
        ylabel('Trotter steps')
    end
    title([num2str(N_m(i)) ' Waveguides']);
    view(-135,30)
    axis tight;
    drawnow;
end

set(h_c_p,'color',0.98*[1 1 1])
drawnow();


A=getframe(h_c_p);
imwrite(A.cdata,'p_convergence_compare.png');
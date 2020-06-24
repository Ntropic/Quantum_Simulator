%Fidelity_Eigenstates_Trotter_Approximation.m
%  -> Determine Eigenstates of exact Unitary
%  -> Fidelity of Trotter approximation for Eigenstates
%  -> Plot approximations
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

%Trotter Iteration number
steps=10;

%Photon numbers Mode numbers
n_max=8;
n_m=[2:2:n_max 3:2:n_max];
N=2; 

%Iteration time
t=pi/4;

%Name
%name=['mwg_trot_compare_fid_n_' num2str(min(n_m)) '_to_' num2str(max(n_m)) '_s_' num2str(steps) '_t_' num2str(t) '_N_' num2str(max(N_m)) '_pauli_ladder'];

connections=[(1:N-1)',(2:N)'];
weights=ones(N-1,1);

figure()
for k=1:length(n_m)
    n=n_m(k);
    %Hamiltonian of the System
    [H]=Gray_Exchange_Hamiltonian_Particles(n); %Mode interaction between 2 modes
    Hij=H2MWG(H,N);                             %Hamiltonian of multiple modes
    %% Exact decomposition
    U_exact=expm(-1i*Hij*t);

    %Subspace of constant particle number n
    indexes=Gray_Indexes_Const_N(n,N); %Subspace dimensions
    U_exact=U_exact(indexes,indexes);

    %Determine the eigenstates and values
    [vec,val]=eig(U_exact);
    val=diag(val);

    for i=1:length(indexes)
        vec_exact{i}=vec(:,i);
    end
    %Now that we have the eigenvalues, we can find the approximations to the
    %eigenvalues of the exact system
    vec_approx=cell(length(indexes),steps);

    %ket_fock_str(indexes,{'fock_gray',n,N})
    %vec

    % Check if vec is symmetric or antisymmetric
    for i=1:length(indexes)
        vec_n=vec(:,i);
        symme=norm(vec_n(1:end)-vec_n(end:-1:1));
        if symme<10^-12
            symmetry(i)=1;
        else
            symmetry(i)=0;
        end
    end

    H_parts=Gray_Hamiltonian_Steps(n); %Gray coding with even uneven terms grouped
    F=[];
    for l=1:steps
        fprintf([' -> ' num2str(l) ' Trotter steps\n'])

        U_approx=U2MWG_Trotter(H_parts,t,l,N); %Create uniteray approximation via Trotter
        U_approx=U_approx(indexes,indexes);    %Subspace

        [min_F(l),F_av(l)]=Min_Av_Fidelity(U_exact,U_approx);
        F(:,l)=Fidelity(U_exact,U_approx,1:length(U_exact),vec)';
    end 

    %% Plotting 
    subplot(2,floor(max(n_m)/2),k)
    fill([1:steps steps:-1:1],[1-min_F,1-F_av(end:-1:1)],[1 1 1]*0.75,'Edgecolor','none');
    hold on;
    for i=1:length(indexes)
        if symmetry(i)==0
            plot(1-F(i,:),'--r')
        else
            plot(1-F(i,:),'--b')
        end
    end
    set(gca,'yscale','log')
    title(['N=' num2str(n)]);
    %xlabel('Trotter steps')
    %ylabel('1-F')
    % o=Order_Of_Convergence(1-F_av_ind(h,i,1:l));
    % plot(o)
    % ylabel('Order of Covergence');
    % xlabel('Trotter steps')
    % title([num2str(N) 'WGs, ' num2str(n) 'Photons'])
    drawnow;
end
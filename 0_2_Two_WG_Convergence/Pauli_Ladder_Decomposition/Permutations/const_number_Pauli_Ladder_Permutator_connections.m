%const_number_Pauli_Ladder_Permutator_connections.m
%Permutates the exchange operations for different total numbers of
%particles
 clear all;
 close all;
 clc;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
clear p;
shortened=fullfile(elem{1:end-3});
addpath(genpath(shortened));

%Photon number
n_max=7;
how_many=100000; %The first how many permutations (maximum)
n_bins=50;
t=pi/4;

mode=0; %0=Ising 1=Boson
norming=0;
for n=3:n_max
    s=ceil(log2(n+1));

    %Create Hamiltonian and seperate into elementary interaction parts
    tree=bintree(s,1);
    gray=all_num_algorithm(tree);

    k=1;
    prefactor=[];
    index=[];
    for i=1:n
        j=n-i+1;
        if mode==0 
            prefactor=[prefactor;1+0*(i*j)];  %Strength of interaction 
        else
            prefactor=[prefactor;sqrt(i*j)];
        end
        gray_num=gray([i+1,j]);
        graybin=dec2bin(gray_num-1,s)';
        graybin=graybin(:)';
        graybin=graybin-'0';

        gray_num2=gray([i,j+1]);
        graybin2=dec2bin(gray_num2-1,s)';
        graybin2=graybin2(:)';
        graybin2=graybin2-'0';
        index(k,:)=[bin2dec(num2str(graybin)) bin2dec(num2str(graybin2))]+1;
        n_i(k)=i-1+j; %Total particles
        m_i(k,:)=[i-1,j,i,j-1]; %Before -> after (each with occupation of both modes)
        k=k+1;
    end

    si=size(index,1);
    %Create elementary Hamiltonians
    H=cell(1,si);
    Ui=H;
    H0=sparse(2^(2*s),2^(2*s));
    Hij=H0;
    for i=1:si
        Hi=H0;
        in=index(i,:);
        Hi(in(1),in(2))=prefactor(i);
        Hi(in(2),in(1))=prefactor(i);
        Hij=Hij+Hi;
        H{i}=Hi;
    end

    Uij=expm(-1i*Hij*t);
    for i=1:si
        Ui{i}=expm(-1i*H{i}*t);
    end

    %Calculate all permutations
    %p=perms(1:si);
    indexes=Gray_Indexes_Const_N(n);
    p=multi_randperm(si,min([factorial(si),how_many]));
    pi=size(p,1);
    F=zeros(1,min([factorial(si),how_many]));
    if min([factorial(si),how_many])>5000
        fprintf([' - > 0/' num2str(min([factorial(si),how_many])) '\n']);
        timer=tic();
        for i=1:min([factorial(si),how_many])
            U=Ui{p(i,1)};
            for j=2:si
                U=U*Ui{p(i,j)};
            end
            [~,F(i)]=Average_Fidelity(U,Uij,indexes);
            if mod(i,500)==0
                fprintf([' - > ' num2str(i) '/' num2str(min([factorial(si),how_many])) ' - t=' num2str(toc(timer)) '\n']);
            end
        end
    else
        for i=1:min([factorial(si),how_many])
            U=Ui{p(i,1)};
            for j=2:si
                U=U*Ui{p(i,j)};
            end
            [~,F(i)]=Average_Fidelity(U,Uij,indexes);
        end
    end

    %%Connections
    diff_p=diff(p,1,2);
    is_conn=abs(diff_p)==1;
    connections=sum(is_conn,2);
    [h c]=hist3([F' connections],'Nbins',[n_bins,n],'CdataMode','auto');
    h_n=[];
    if norming==1
        for i=1:size(h,2)
            if sum(h(:,i))>0
                h_n(:,i)=h(:,i)/sum(h(:,i));
            else
                h_n(:,i)=h(:,i);
            end
        end
    else
        h_n=h;
    end
    d=c{1};
    figure()
    imagesc(0:n-1,d(end:-1:1),h_n)%(end:-1:1,:));
    set(gca,'XTick',0:n-1);
    xlabel('\# of Connections')
    ylabel('$F_1$')
    title(['n=' num2str(n)])
    drawnow()
    matlab2tikz(['Connections\const_n_connections_' num2str(n) '.tex'],'standalone',true,'parseStrings',false,'width','6.2cm','height','4.5cm');

%     n_c=[];
%     for i=0:n-1
%         n_c(i+1)=sum(connections==i);
%     end
%     figure()
%     plot(0:n-1,n_c)
end
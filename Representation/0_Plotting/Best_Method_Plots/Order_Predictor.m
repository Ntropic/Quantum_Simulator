function [ F_pred ] = Order_Predictor( F )
%ORDER_PREDICTOR predicts order of Convergence for Fidelities 
% F_pred declares first the order of convergence and 2nd the error at
% maximum steps (last dimension) and 3rd the maximum step number

ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0.8 2];
    
if isa(F,'cell')
    F_pred=cell(3,length(F));
    for i=1:length(F)
        F_i=F{i};
        sizes=size(F_i);
        steps=sizes(end);
        dim=length(sizes);
        F_i=permute(F_i,[dim (1:(dim-1))]); %Put Trotter steps to the beginning
        
        %how_many fits?
        how_many=prod(size(F_i))/steps;
        indexes=1:steps:prod(size(F_i))+1;
        F_i=F_i(:);
        F_pred_order=zeros(how_many,1);
        F_pred_fn=zeros(how_many,1);
        F_pred_n=zeros(how_many,1);
        for j=1:how_many
            f=1-F_i(indexes(j):(indexes(j+1)-1)); %error rates for current fit
            o=Order_Of_Convergence(f,(1:steps)+0.5);
            [fitresult, gof] = fit( (1:steps-1)', o', ft, opts );
            fitc=fitresult.c;   %Most relevant fit parameter
            
            if fitc<0
                opts.StartPoint = [0 0.8 o(end)];
                [fitresult, gof] = fit( (1:steps-1)', o', ft, opts );
                fitc=fitresult.c;
                if fitc<0
                    fitc=o(end);
                end
            end
                
            F_pred_order(j)=fitc;
            F_pred_fn(j)=f(steps-1);
            F_pred_n(j)=steps-1;
        end
        F_pred_order=reshape(F_pred_order,sizes(1:end-1));
        F_pred_fn=reshape(F_pred_fn,sizes(1:end-1));
        F_pred_n=reshape(F_pred_n,sizes(1:end-1));
        F_pred{1,i}=F_pred_order;
        F_pred{2,i}=F_pred_fn;
        F_pred{3,i}=F_pred_n;
    end
else
    sizes=size(F);
    steps=sizes(end);
    dim=length(sizes);
    F=permute(F,[dim (1:(dim-1))]); %Put Trotter steps to the beginning

    %how_many fits?
    how_many=prod(size(F))/steps;
    indexes=1:steps:prod(size(F_i));
    F=F(:);
    F_pred_order=zeros(how_many,1);
    F_pred_fn=zeros(how_many,1);
    F_pred_n=zeros(how_many,1);
    for j=1:how_many
        f=1-F(indexes(j):(indexes(j+1)-1)); %error rates for current fit
        o=Order_Of_Convergence(f,1:steps);
        [fitresult, gof] = fit( 1:steps-1, o, ft, opts );
        fitc=fitresult.c;   %Most relevant fit parameter

        F_pred_order(j)=fitc;
        F_pred_fn(j)=f(steps-1);
        F_pred_n(j)=steps-1;
    end
    F_pred_order=reshape(F_pred_order,sizes(1:end-1));
    F_pred_fn=reshape(F_pred_fn,sizes(1:end-1));
    F_pred_n=reshape(F_pred_n,sizes(1:end-1));
    F_pred{1}=F_pred_order;
    F_pred{2}=F_pred_fn;
    F_pred{3}=F_pred_n;
end
end


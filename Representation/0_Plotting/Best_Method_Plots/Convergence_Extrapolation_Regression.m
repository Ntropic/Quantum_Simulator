function [fitresult, gof] = Convergence_Extrapolation_Regression(o,n,plotting)
%Extrapolate the Convergence Order from Existing Data with the Regression
%model O(n)=a*exp(-b*n)+c -> converges to c for positive b
    if nargin==1
        [xData, yData] = prepareCurveData( [], o );
    else
        [xData, yData] = prepareCurveData( n, o );
    end
    if nargin<3
        plotting=0;
    end

    % Set up fittype and options.
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [-1 0.8 2];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    if plotting==1
        % Plot fit with data.
        figure( 'Name','Order Of Convergence');
        h = plot( fitresult, xData, yData );
        legend( h, 'o', '\mathcal{O}(n)', 'Location', 'NorthEast' );
        % Label axes
        ylabel('\mathcal{O}');
        xlabel('n');
        grid on
    end
end


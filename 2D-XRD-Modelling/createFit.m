function [fitresult, gof] = createFit(t, f)
%CREATEFIT(T,F)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : t
%      Y Output: f
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t, f );

% Set up fittype and options.
ft = 'pchipinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
%h = plot( fitresult, xData, yData );
%legend( h, 'f vs. t', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
%xlabel( 't' );
%ylabel( 'f' );
grid on



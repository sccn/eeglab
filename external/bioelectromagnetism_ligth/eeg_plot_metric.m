function eeg_plot_metric

% eeg_plot_metric - create TeX plot label (Y is uV, X is msec)

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2000, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Font.FontName   = 'Helvetica';
Font.FontUnits  = 'Pixels';
Font.FontSize   = 12;
Font.FontWeight = 'normal';
Font.FontAngle  = 'normal';

ylab = get(gca,'ylabel');
set(ylab,'string','\muV','Rotation',0,Font,'units','normalized',...
    'Position',[-0.10  0.95 0]);
xlab = get(gca,'xlabel');
set(xlab,'string','msec','Rotation',0,Font,'units','normalized',...
    'Position',[ 1.07 -0.01 0]);
    
return

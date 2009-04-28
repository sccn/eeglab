function[p] = gui_erp_plot(p)

% gui_erp_plot - eeg_toolbox plot for ERP data
%
%[p] = erp_plot(p)
%
% creates a GUI interface for ERP measurement and
% topographic mapping functions.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:56 $

% Licence:  GNU GPL, no express or implied warranties
% History:  10/2003, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plotfig = figure('Name',p.volt.file,...
  'NumberTitle','off',...
  'UserData',p);

movegui(plotfig,'center');

plot(p.volt.timeArray,p.volt.data);
axis tight;
eeg_plot_metric;

[Xpoint,Ypoint] = eeg_crosshair('init',p);

return

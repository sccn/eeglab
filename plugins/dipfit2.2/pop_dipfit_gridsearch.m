% pop_dipfit_gridsearch() - scan all ICA components with a single dipole
%     on a regular grid spanning the whole brain. Any dipoles that explains
%     a component with a too large relative residual variance is removed.
%
% Usage: 
%  >> EEGOUT = pop_dipfit_gridsearch( EEGIN ); % pop up interactive window
%  >> EEGOUT = pop_dipfit_gridsearch( EEGIN, comps );
%  >> EEGOUT = pop_dipfit_gridsearch( EEGIN, comps, xgrid, ygrid, zgrid, thresh )
%
% Inputs:
%   EEGIN     - input dataset
%   comps     - [integer array] component indices
%   xgrid     - [float array] x-grid. Default is 10 elements between
%               -1 and 1.
%   ygrid     - [float array] y-grid. Default is 10 elements between
%               -1 and 1.
%   zgrid     - [float array] z-grid. Default is 10 elements between
%               -1 and 1.
%   thresh    - [float] threshold in percent. Default 40.
%
% Outputs:
%   EEGOUT      output dataset
%
% Authors: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%          Arnaud Delorme, SCCN, La Jolla 2003
%          Thanks to Nicolas Robitaille for his help on the CTF MEG
%          implementation

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl/

% Copyright (C) 2003 Robert Oostenveld, SMI/FCDC roberto@smi.auc.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EEGOUT, com] = pop_dipfit_gridsearch(EEG, select, xgrid, ygrid, zgrid, reject );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
  help pop_dipfit_gridsearch;
  return;
end;

EEGOUT = EEG;
com = '';

if ~isfield(EEG, 'chanlocs')
  error('No electrodes present');
end

if ~isfield(EEG, 'icawinv')
  error('No ICA components to fit');
end

if ~isfield(EEG, 'dipfit')
  error('General dipolefit settings not specified');
end

if ~isfield(EEG.dipfit, 'vol') & ~isfield(EEG.dipfit, 'hdmfile')
  error('Dipolefit volume conductor model not specified');
end

dipfitdefs
if strcmpi(EEG.dipfit.coordformat, 'CTF')
    maxrad = 8.5;
    xgridstr     = sprintf('linspace(-%2.1f,%2.1f,11)', maxrad, maxrad);
    ygridstr     = sprintf('linspace(-%2.1f,%2.1f,11)', maxrad, maxrad);
    zgridstr     = sprintf('linspace(0,%2.1f,6)', maxrad);
end;
if nargin < 2
  % get the default values and filenames
  promptstr = { 'Component(s) (not faster if few comp.)', ...
      'Grid in X-direction', ...
      'Grid in Y-direction', ...
      'Grid in Z-direction', ...
      'Rejection threshold RV(%)' };
  
  inistr    = { 
    [ '1:' int2str(size(EEG.icawinv,2)) ], ...
      xgridstr, ...
      ygridstr, ...
      zgridstr, ...
      rejectstr };
    
    result = inputdlg2( promptstr, 'Batch dipole fit -- pop_dipfit_gridsearch()', 1,  inistr, 'pop_dipfit_gridsearch');
    
    if length(result)==0
      % user pressed cancel
      return
    end
    
    select = eval( [ '[' result{1} ']' ]);
    xgrid  = eval( result{2} );
    ygrid  = eval( result{3} );
    zgrid  = eval( result{4} );
    reject = eval( result{5} ) / 100;	% string is in percent
    options = { };
  else
    if nargin < 2
      select = [1:size(EEG.icawinv,2)];
    end;
    if nargin < 3
      xgrid  = eval( xgridstr );
    end;
    if nargin < 4
      ygrid  = eval( ygridstr );
    end;
    if nargin < 5
      zgrid  = eval( zgridstr );
    end;
    if nargin < 6
      reject  = eval( rejectstr );
    end;
    options = { 'waitbar' 'none' };
  end;
  
  % perform batch fit with single dipole for all selected channels and components
  % warning off;
  warning backtrace off;  
  EEGOUT = dipfit_gridsearch(EEG, 'component', select, 'xgrid', xgrid, 'ygrid', ygrid, 'zgrid', zgrid, options{:});
  warning backtrace on;
  EEGOUT.dipfit.model  = dipfit_reject(EEGOUT.dipfit.model, reject);

  % FIXME reject is not being used at the moment
  disp('Done');
  com = sprintf('%s = pop_dipfit_gridsearch(%s, %s);', ...
                inputname(1), inputname(1), vararg2str( { select xgrid, ygrid, zgrid reject }));
  

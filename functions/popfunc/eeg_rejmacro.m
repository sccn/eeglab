% eeg_rejmacro() - Internal EEGLAB macro for all pop_ functions that
%                  perform data rejection.   
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 03-08-02 include past rejections in eegplot -ad

% this macro is used by pop_eegplot, pop_eegthresh, pop_const, pop_jointprob, pop_eegthreshfreq, pop_kurtosis

% inputs:
%   rej - input rejection in this workspace (string)
%   rejE - electrode rejection in this workspace (string)
%   elecrange - range of selected electrodes
%   reject - superpose with previous rejections
% optional : nbpnts = number of points per trial

% outputs:
%   rejection arrays ('oldrej', 'oldrejE'), 'command' for eeglot

if ~exist('elecrange')
	help eeg_rejmacro;
	error('Error: eeg_rejmacro cannot be called from the command line');
end;	

if ~exist('nbpnts') nbpnts = EEG.pnts; end;

% mix all type of rejections
% --------------------------
if superpose == 2
    com2 = ['if ~isempty(TMPREJ), ' ...
            '  icaprefix = ' fastif(icacomp, '''''', '''ica''') ';' ...
            '  for indextmp = 1:length(EEG.reject.disprej),' ...
            '     eval([ ''colortmp = EEG.reject.rej'' EEG.reject.disprej{indextmp} ''col;'' ]);' ... 
            '     [tmprej tmprejE] = eegplot2trial(TMPREJ,' int2str(nbpnts) ', EEG.trials, colortmp, []);' ...
            '     if ~isempty(tmprejE),' ...
            '          tmprejE2 = zeros(EEG.nbchan, length(tmprej));' ...
            '          tmprejE2([' int2str(elecrange) '],:) = tmprejE;' ... 
            '     else,' ...
            '          tmprejE2 = [];' ...
            '     end;' ...
            '     eval([ ''EEG.reject.'' icaprefix ''rej'' EEG.reject.disprej{indextmp} ''= tmprej;'' ]);' ...
            '     eval([ ''EEG.reject.'' icaprefix ''rej'' EEG.reject.disprej{indextmp} ''E = tmprejE2;'' ]);' ...
            '  end;' ];
%           '     colortmp, tmpstr = [ ''EEG.reject.'' icaprefix ''rej'' EEG.reject.disprej{indextmp} ], sum(tmprej),' ...
else      
    com2 = [ 'if ~isempty(TMPREJ),' ...
             '  icaprefix = ' fastif(icacomp, '''''', '''ica''') ';' ...
             '  [tmprej tmprejE] = eegplot2trial(TMPREJ,' int2str(nbpnts) ', EEG.trials, [' num2str(colrej) '], []);' ...
             '  if ~isempty(tmprejE),' ...
             '     tmprejE2 = zeros(EEG.nbchan, length(tmprej));' ...
             '     tmprejE2([' int2str(elecrange) '],:) = tmprejE;' ... 
             '  else,' ...
             '     tmprejE2 = [];' ...
             '  end;' ...
             macrorej '= tmprej;' macrorejE '= tmprejE2;' ...
             ... % below are manual rejections
             '  tmpstr = [ ''EEG.reject.'' icaprefix ''rejmanual'' ];' ...
             '  if ~isempty(tmprej) eval([ ''if ~isempty('' tmpstr ''),'' tmpstr ''='' tmpstr ''| tmprej; else, '' tmpstr ''=tmprej; end;'' ]); end;' ...
             '  if ~isempty(tmprejE2) eval([ ''if ~isempty('' tmpstr ''E),'' tmpstr ''E='' tmpstr ''E| tmprejE2; else, '' tmpstr ''E=tmprejE2; end;'' ]); end;' ];
%             '  size(tmprejE2), eval([''disp(size('' tmpstr ''E))'']),' ...
end;
% text commented below to fix BUG 478
%             '  [tmprej tmprejE] = eegplot2trial(TMPREJ,' int2str(nbpnts) ', EEG.trials, [' num2str(EEG.reject.rejmanualcol) '], []);' ...
%             '  if ~isempty(tmprejE),' ...
%             '     tmprejE2 = zeros(EEG.nbchan, length(tmprej));' ...
%             '     tmprejE2([' int2str(elecrange) '],:) = tmprejE;' ... 
%             '  else,' ...
%             '     tmprejE2 = [];' ...
%             '  end;' ...
if reject
    com2 = [ com2 ...
           '[EEGTMP LASTCOM] = pop_rejepoch(EEG, tmprej, 1);' ...
		   'if ~isempty(LASTCOM),'...
           '   [ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
           '   if ~isempty(tmpcom),' ... 
           '     EEG = eegh(LASTCOM, EEG);' ...
           '     eegh(tmpcom);' ...
           '     eeglab(''redraw'');' ...
            '  end;' ...
	       'end; clear EEGTMP tmpcom; end;' ] ;
else
	com2 = [com2 ...
			' warndlg2(strvcat(''Epochs (=trials) marked for rejection have been noted.'',' ...
			'''To actually reject these epochs, use '', ''Tools > Reject data epochs > Reject marked epochs''), ''Warning'');' ...
			'[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eeglab(''redraw''); end;' ];
end; 
if ~exist('topcommand')
	topcommand = [];
end;

% the first part is used to convert the eegplot output
command = [  com2 topcommand 'clear indextmp colortmp icaprefix tmpcom tmprej tmprejE tmprejE2 TMPREJ;' ];

if all(colrej == EEG.reject.rejmanualcol)
	oldrej = [];  % for manual rejection, old rejection are
	oldrejE = []; % the current rejection
else
	oldrej  = eval(macrorej);
	oldrejE = eval(macrorejE);
end;

switch superpose
 case 0, rejeegplot = trial2eegplot(  rej, rejE, nbpnts, colrej);
 case 1, rejeegplottmp = trial2eegplot(  oldrej, oldrejE, nbpnts, min(colrej+0.15, [1 1 1]));
         if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplottmp ]; 
		 else rejeegplot = []; end;
         rejeegplottmp = trial2eegplot(  rej, rejE, nbpnts, colrej);
         if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplot; rejeegplottmp ]; end;
 case 2, 
  rejeegplot = [];
  for index = 1:length(EEG.reject.disprej)
	  if ~isempty(EEG.reject.disprej{index})
		  eval([ 'colortmp = EEG.reject.rej' EEG.reject.disprej{index} 'col;']); 
		  if any(colortmp ~= colrej) % test if current rejection (if color different)
			  if icacomp == 0 % ica
				  currentname = [ 'EEG.reject.icarej' EEG.reject.disprej{index} ];
			  else
				  currentname = [ 'EEG.reject.rej' EEG.reject.disprej{index} ];
			  end;
			  currentcolor =  [ 'EEG.reject.rej' EEG.reject.disprej{index} 'col' ];
			  %if strcmp(EEG.reject.disprej{index}, 'manual')
			  %	  currentcolor = [ 'min(' currentcolor '+0.15, [1 1 1])' ];
			  %end; % using this test, manual rejections won't be added to current rej
			  eval( [ 'rejeegplottmp = trial2eegplot( ' currentname ',' currentname ...
					  'E, nbpnts,' currentcolor ');' ]);
			  if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplot; rejeegplottmp ]; end;
		  end;
	  end;
  end;
  rejeegplottmp = trial2eegplot(  rej, rejE, nbpnts, colrej);
  if ~isempty(rejeegplottmp), rejeegplot = [ rejeegplot; rejeegplottmp ]; end;
end;
if ~isempty(rejeegplot)
	rejeegplot = rejeegplot(:,[1:5,elecrange+5]);
else
	rejeegplot = [];
end;
eegplotoptions = { 'events', EEG.event, 'winlength', 5, 'winrej', ...
				   rejeegplot, 'xgrid', 'off', 'wincolor', EEG.reject.rejmanualcol, ...
				   'colmodif', { { EEG.reject.rejmanualcol EEG.reject.rejthreshcol EEG.reject.rejconstcol ...
                                   EEG.reject.rejjpcol     EEG.reject.rejkurtcol   EEG.reject.rejfreqcol } } };

if ~isempty(EEG.chanlocs) & icacomp == 1
    if exist('elecrange')
        eegplotoptions = { eegplotoptions{:}  'eloc_file', EEG.chanlocs(elecrange) };
    else
        eegplotoptions = { eegplotoptions{:}  'eloc_file', EEG.chanlocs };
    end;
else 
    if exist('elecrange')
        for index = 1:length(elecrange)
            tmpstruct(index).labels = int2str(elecrange(index));
        end;
        eegplotoptions = { eegplotoptions{:}  'eloc_file' tmpstruct };
    end;
end;
if ~reject
	eegplotoptions = { eegplotoptions{:}  'butlabel', 'UPDATE MARKS' };
end;

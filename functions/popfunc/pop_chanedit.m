% pop_chanedit() - Edit channel locations.
%
% Usage: >> newchans = pop_chanedit( chans, 'key1', value1, ...
%                                           'key2', value2, ... );
% Input:
%   chans - channel EEGLAB structure
%
% Optional inputs:
%   'convert'     - { conversion_type args } conversion type is 'cart2topo'
%                   'sph2topo', 'topo2sph', 'sph2cart', 'cart2sph', 'chancenter'. 
%                   Help can be found in the function having the same name. Args are
%                   only relevant for 'chancenter'.
%   'transform'     - string command for manipulating arrays. 'chan' is a full 
%                   electrode info. Fields can be manipulated using 'labels', 'theta'
%                   'radius' (polar angle and radius), 'X', 'Y', 'Z' (cartesian 3D) or
%                   'sph_theta', 'sph_phi', 'sph_radius' for spherical horizontal angle, 
%                   azimut and radius. Ex: 'chans(3) = chans(14)'. 'X = -X' or a more complex 
%                   transform separated by ';': 'TMP = X; X = Y; Y = TMP'
%   'changechan'  - {num value1 value2 value3 ...} Change the values of
%                   all fields the channel.
%   'changefield' - {num field value} change field value for channel number num. 
%                   (Ex: {34 'theta' 320.4}).
%   'add'         - {num label theta radius X Y Z sph_theta sph_phi sph_radius } 
%                   Insert channel before channel number num with the specified values. 
%                   If the number of values if less than ten, fields are filled with 0.
%   'delete'      - vector of indices of channel to delete
%   'shrink'      - topographical polar shrink factor (see >> help topoplot) 
%   'load'        - [filename|{filename, 'key', 'val'}] load channel location file
%                   optional arguments (such as file format) to the function readlocs()
%                   can be given if the input is a cell array.
%   'save'        - 'filename' save text file with channel info
%
% Outputs:
%   newchans      - channel EEGLAB structure
%
% Ex:  EEG = pop_chanedit(EEG,'load', { 'dummy.elp' 'elp' }, 'delete', [3 4], ...
%          'convert', { 'xyz->polar' [] -1 1 }, 'save', 'myfile.loc' )
%        % load polhemus file, delete two channels, convert to polar (see 
%        % cart2topo() for arguments) and save into 'myfile.loc'.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 20 April 2002
%
% See also: readlocs()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 15 March 2002, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.28  2002/10/15 17:06:44  arno
% drawnow
%
% Revision 1.27  2002/10/15 17:02:38  arno
% drawnow
%
% Revision 1.26  2002/10/10 21:17:30  arno
% debug command line call and display (when deleting last channel)
%
% Revision 1.25  2002/09/23 17:57:25  arno
% debug shrink off
%
% Revision 1.24  2002/08/12 21:38:24  arno
% text
%
% Revision 1.23  2002/08/12 21:17:59  arno
% debug
%
% Revision 1.22  2002/08/12 18:50:22  arno
% errordlg2
%
% Revision 1.21  2002/08/12 18:30:05  arno
% quesdlg2
%
% Revision 1.20  2002/08/12 16:43:01  arno
% debug inputdlg2
%
% Revision 1.19  2002/08/12 16:37:24  arno
% inputdlg2
%
% Revision 1.18  2002/08/02 15:41:23  arno
% adding file format input
%
% Revision 1.17  2002/07/30 00:40:48  arno
% debugging shrink
%
% Revision 1.16  2002/06/25 14:27:40  arno
% typo and 3d plot check
%
% Revision 1.15  2002/05/21 20:44:31  scott
% removed ; from evalin() calls -sm
%
% Revision 1.14  2002/05/03 00:49:24  arno
% debugging channel center
%
% Revision 1.13  2002/05/02 23:37:27  scott
% editting -sm
%
% Revision 1.12  2002/05/02 23:35:15  scott
% formula -> transform -sm & ad
%
% Revision 1.11  2002/05/02 23:29:56  scott
% Num -> number -sm
%
% Revision 1.10  2002/05/02 23:18:27  scott
% editting -sm
%
% Revision 1.9  2002/05/02 23:16:44  scott
% editting -sm
%
% Revision 1.8  2002/05/02 23:12:17  arno
% editing text
%
% Revision 1.7  2002/05/02 01:46:03  arno
% debugging for full consistency
%
% Revision 1.6  2002/05/01 19:38:38  arno
% returning shrink factor
%
% Revision 1.5  2002/05/01 03:30:55  arno
% editing interface
%
% Revision 1.4  2002/05/01 03:29:28  arno
% correct typo
%
% Revision 1.3  2002/05/01 03:24:27  arno
% function checkchans
%
% Revision 1.2  2002/05/01 03:22:22  arno
% testelp->plotchans3d
%
% Revision 1.1  2002/05/01 03:14:50  arno
% Initial revision
%

% hidden parameter
%   'gui' - [figure value], allow to process the same dialog box several times

function [chans, com] = pop_chanedit(chans, varargin);

com ='';
if nargin < 1
   help pop_editeventvals;
   return;
end;	

nbchan = length(chans);
[chans shrinkfact]= checkchans(chans);

allfields = fieldnames(chans);
if nargin < 2
	totaluserdat = {};
	ingui = 1;
	guimodif = 0;
	while ingui
		commentfields = { 'Channel label ("label")', 'Polar angle ("theta")', 'Polar radius ("radius")', 'Cartesian X ("X")', 'Cartesian Y ("Y")', 'Cartesian Z ("Z")', ...
						  'Spherical horiz. angle ("sph_theta")', 'Spherical azimuth angle ("sph_phi")', 'Spherical radius ("sph_radius")' };
		% transfer channel to global workspace
		global chantmp;
		chantmp = chans;
		evalin('base', [ 'global chantmp ' ]);
		
		% add field values
		% ----------------
		geometry = { 1 };
		tmpstr = sprintf('Channel information ("field_name"):');
		uilist = { { 'Style', 'text', 'string', tmpstr, 'fontweight', 'bold'  } };
		endgui = 'set(findobj(''parent'', gcbf, ''tag'', ''ok''), ''userdata'', ''stop'');';
		transform = ['inputdlg2({''Enter transform: (Ex: TMP=X; X=-Y; Y=TMP or Y(3) = X(2), etc.'' }, ''Transform'', 1, { '''' }, ''pop_chanedit'');'];
		uiconvert = { { 'Style', 'pushbutton', 'string', '3D center', 'callback', ...
						['comtmp = {''convert'' {''chancenter'' [] 1}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'xyz->polar', 'callback', ...
						['comtmp = {''convert'' {''cart2topo'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'sph->polar', 'callback', ...
						['comtmp = {''convert'' {''sph2topo'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'polar->sph', 'callback', ...
						['comtmp = {''convert'' {''topo2sph'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'sph->xyz'  , 'callback', ...
						['comtmp = {''convert'' {''sph2cart'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'xyz->sph'  , 'callback', ...
						['comtmp = {''convert'' {''cart2sph'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'Transform axes', 'callback', ...
						['comtmp = {''transform'' ' transform '};' endgui] }, ...
					  {} { } { } };
		%{ 'Style', 'pushbutton', 'string', 'UNDO LAST ', 'callback', '' } { } { } };
		for index = 1:length(allfields)
			geometry = { geometry{:} [1.5 1 0.2 1] };
			uilist   = { uilist{:}, ...
						 { 'Style', 'text', 'string', commentfields{index} }, ...
						 { 'Style', 'edit', 'tag', [ 'chanedit' allfields{index}], 'string', num2str(getfield(chans,{1}, allfields{index})), 'callback', ...
						   [ 'valnum   = str2num(get(findobj(''parent'', gcbf, ''tag'', ''chaneditnumval''), ''string''));' ...
							 'editval = get(gcbo, ''string'');' ...
							 'if ~isempty(str2num(editval)), editval = str2num(editval);  end;' ...
							 'eval([ ''chantmp(valnum).' allfields{index} '= editval;'']);' ...
							 'olduserdat = get( gcbf, ''userdata'');' ...
							 'if isempty(olduserdat), olduserdat = {}; end;' ...
							 'set( gcbf, ''userdata'', { olduserdat{:} ''changefield'' { valnum ''' allfields{index} ''' editval }});' ...
							 'clear editval valnum olduserdat;' ] } { } uiconvert{index} };
		end;
		
		% add buttons
		% -----------
		geometry =  { geometry{:} [1] [1 0.7 0.4 1.8 0.2 0.7 1] [1 0.7 0.7 1 0.7 0.7 1] };
		callpart1 = [ 'valnum   = str2num(get(findobj(''tag'', ''chaneditnumval''), ''string''));' ];
		callpart2 = [ 'set(findobj(''tag'', ''chaneditnumval''), ''string'', num2str(valnum));' ];
		for index = 1:length(allfields)
			callpart2 = [ callpart2  'set(findobj(''tag'', ''chanedit' allfields{index} ...
						  '''), ''string'', num2str(chantmp(valnum).' allfields{index} '));' ];
		end;
		callpart2 = [ callpart2 'set(findobj(''tag'', ''chaneditscantitle''), ' ...
					  '''string'', [''Channel number (of '' int2str(length(chantmp)) '')'']);' ...
                      'set(findobj(''tag'', ''chaneditnumval''), ''string'', int2str(valnum));'  ]; 
		callpart2 = [ callpart2 'clear orivalnum valnum;' ];
		
		uilist   = { uilist{:}, ...
					 { }, ...
					 { },{ },{ }, {'Style', 'text', 'string', ['Channel number (of ' int2str(length(chans)) ')'], ...
					'fontweight', 'bold', 'tag', 'chaneditscantitle' }, { },{ },{ }, ...
					 { 'Style', 'pushbutton', 'string', 'Delete chan',  'callback', [callpart1 'orivalnum = valnum; chantmp(valnum) = [];' ...
					'valnum = min(valnum,length(chantmp));' ...
                    'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
					'set( gcbf, ''userdata'', { olduserdat{:} ''delete'', orivalnum }); clear olduserdat' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '<<', 'callback', [callpart1 'valnum = max(valnum-10,1);' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '<',  'callback', [callpart1 'valnum = max(valnum-1,1);' callpart2 ] }, ...
					 { 'Style', 'edit', 'string', '1', 'tag', 'chaneditnumval', 'callback', ...
					   [callpart1 'valnum = min(str2num(get(gcbo, ''string'')),length(chantmp));' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '>',  'callback', [callpart1 'valnum = min(valnum+1,length(chantmp));' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '>>', 'callback', [callpart1 'valnum = min(valnum+10,length(chantmp));' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', 'Insert chan',  'callback', [callpart1 ...
					'chantmp(end+3) = chantmp(end);' ...
					'chantmp(valnum+1:end-2) = chantmp(valnum:end-3);' ...
					'chantmp(valnum) = chantmp(end-1);' ...
					'if isfield(chantmp, ''epoch''), chantmp(valnum).epoch = chantmp(valnum+1).epoch; end;' ...
					'chantmp = chantmp(1:end-2);' ...
					'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
					'tmpcell = cell(1,1+length(fieldnames(chantmp))); tmpcell{1} =valnum;' ...
					'set( gcbf, ''userdata'', { olduserdat{:} ''add'', tmpcell }); clear tmpcell olduserdat' callpart2 ] }, ...
				   };
		
		% add sorting options
		% -------------------
		if isstr(shrinkfact)
			evalin('base', ['tmpshrink = ''' shrinkfact ''';' ]);
		else
			evalin('base', ['tmpshrink = ' num2str(shrinkfact) ';' ]);
		end;
		plot2dcom = [ 'tmpshrink = ''off'';' ...
					  'if get(findobj(''parent'', gcbf,  ''string'', ''Auto shrink''), ''value'')' ...
					  '   tmpshrink = ''force'';' ...
					  'else ' ...
					  '   if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''shrinkfactor''), ''string''))' ...
					  '       tmpshrink = str2num(get(findobj(''parent'', gcbf, ''tag'', ''shrinkfactor''), ''string''));' ...
					  '   end;' ...
					  'end;' ...
					  'figure; topoplot([],chantmp, ''style'', ''blank'', ''electrodes'', ''labelpoint'', ''shrink'', tmpshrink);' ];
		geometry = { geometry{:} [1] [1 1.4 0.6 1 1] [1] [1 1 1 1]};
		uilist = {  uilist{:},...
					{ } ...
					{ 'Style', 'pushbutton', 'string', 'Plot 2D', 'callback', plot2dcom },... 
					{ 'Style', 'text', 'string', '2D shink factor (-1 to 1)'} ...
					{ 'Style', 'edit', 'string', fastif(isnumeric(shrinkfact), num2str(shrinkfact), ''), ...
					  'tag', 'shrinkfactor', 'callback', 'tmpshrink=get(gcbo, ''string''); set(findobj(''tag'', ''shrinkbut''), ''value'', 0);' } ...
					{ 'Style', 'checkbox', 'tag', 'shrinkbut', 'string', 'Auto shrink', 'value', strcmp(shrinkfact, 'force'), ...
					  'callback', 'if get(gcbo, ''value''), tmpshrink=''force''; else tmpshrink=''off''; end; set(findobj(''tag'', ''shrinkfactor''), ''string'', '''');' } ...
					{ 'Style', 'pushbutton', 'string', 'Plot 3D (XYZ)', 'callback', [ 'if ~isempty(chantmp(1).X),' ...
					 'plotchans3d([cell2mat({chantmp.X})'' cell2mat({chantmp.Y})'' cell2mat({chantmp.Z})''],' ...
					'{chantmp.labels}); else disp(''cannot plot: no XYZ coordinates''); end;'] } ...
					{}, { 'Style', 'pushbutton', 'string', 'Read locations', 'callback', ...
						  ['[tmpf tmpp] = uigetfile(''*'', ''Load a channel location file''); drawnow;' ...
						   'tmpfmt = inputdlg2({[''File format [loc|sph|xyz|polhemus|besa|chanedit] ([]=use extension)'' ]}, ''File format'', 1, { '''' }, ''readlocs'');' ...
						   'if isempty(tmpfmt), tmpfmt = {[]}; tmpp = []; tmpf = []; end;' ...
						   'comtmp = {''load'' { [tmpp tmpf ] ''filetype'' tmpfmt{1} } }; clear tmpfmt;' endgui ] }, ...
					{ 'Style', 'pushbutton', 'string', 'Read help', 'callback', 'pophelp(''readlocs.m'');' }, ...	
					{ 'Style', 'pushbutton', 'string', 'Save ASCII', 'callback', ...	
					  ['[tmpf tmpp] = uiputfile(''*'', ''Save channel location file''); drawnow' ...
					   'comtmp = {''save'' [tmpp tmpf] };' endgui ] }, ...
					{ 'Style', 'pushbutton', 'enable', 'off', 'string', 'Save BESA' }, ...	
				 };
		
		if guimodif
			set(currentfig, 'userdat', {});
			eval([ callpart1 callpart2 ]);
			[results userdat returnmode] = inputgui( geometry, uilist, 'pophelp(''pop_chanedit'');', 'Edit channel info -- pop_chanedit()', {}, currentfig );
		else 
			[results userdat returnmode] = inputgui( geometry, uilist, 'pophelp(''pop_chanedit'');', 'Edit channel info -- pop_chanedit()', {}, 'noclose' );
			currentfig = gcf;
		end; 
		if length(results) == 0, return; end;
		
		% transfer events back from global workspace
		chans = evalin('base', 'chantmp');
		%evalin('base', 'clear chantmp');
		totaluserdat = { totaluserdat{:} userdat{:}};
		if evalin('base', 'exist(''comtmp'')') == 1
			tmpcom = evalin('base', 'comtmp');
			try, 
				chans = pop_chanedit(chans, tmpcom{:}); % apply modification to channel structure
                if iscell(tmpcom{2}) & (length(tmpcom{2}) == 2) & isstr(tmpcom{2}{2}) & strcmpi(tmpcom{2}{2}, 'gui'),
                    tmpcom = { tmpcom{1} tmpcom{2}{1} };
                end;
                totaluserdat = { totaluserdat{:} tmpcom{:} };
			catch
				errordlg2(lasterr, 'Channel location error');
				returnmode = 'no';
			end;	
			evalin('base', 'clear comtmp');
			chans = checkchans(chans);
		end;
		
		% handle arguments
		% ----------------
		if strcmp(returnmode, 'retuninginputui')
			ingui = 0;
			if nbchan ~= 0 & nbchan ~= length(chans)
				if popask(strvcat(['The number of channel (' int2str(length(chans)) ') does not correspond to the'], ...
								  ['initial number of channel (' int2str(nbchan) '), so the channel information'], ...
								  'may be removed if this function was called from EEGLAB'))
					ingui = 1;
				end;
			end;	
		else 
			guimodif = 1;
		end;
	end;
	
	% not in the loop, returning
	% --------------------------
	if ~isempty(findobj('parent', gcf, 'tag','shrinkfactor'))
		if evalin('base', 'exist(''tmpshrink'')') == 1
			tmpshrink = evalin('base', 'tmpshrink');
			evalin('base', 'clear tmpshrink');
			if ~strcmp(num2str(tmpshrink), num2str(shrinkfact))
				if ~isempty(str2num(tmpshrink))  chans(1).shrink = str2num(tmpshrink);
				else                             chans(1).shrink = tmpshrink;
				end;
				totaluserdat{end+1} = 'shrink';
				totaluserdat{end+1} = chans(1).shrink;
			end; 
		end;
		close(gcf);
		if ~isempty( totaluserdat )
			if isempty(inputname(1)), varname = 'EEG.chanlocs';
			else varname = inputname(1);
			end;
            com = sprintf('%s=pop_chanedit(%s, %s);', varname, varname, vararg2str(totaluserdat));
		end;
	end;
	evalin('base', 'clear global tmpchan; clear tmpchan;');
else 
	 args = varargin;
	 % no interactive inputs
	 % scan all the fields of g
	 % ------------------------
	 for curfield = 1:2:length(args)
		 switch lower(args{curfield})
		  case 'convert', 
		   if iscell(args{curfield+1})
			   method=args{curfield+1}{1};
			   extraargs = args{curfield+1}(2:end);
		   else
			   method=args{curfield+1};
			   extraargs = {''};
		   end;
		   if isstr(extraargs{1}) & strcmp(extraargs{1}, 'gui') & ~strcmp(method, 'chancenter')
			   tmpButtonName=questdlg2( strvcat('This will modify fields in the channel structure', ...
					'Are you sure you want to apply this function ?'), 'Confirmation', 'Cancel', 'Yes','Yes');
			   if ~strcmp(tmpButtonName, 'Yes'), return; end;
		   end;
		   switch method
			case 'chancenter',
			 extraargs
			 if isempty(extraargs)
				 [X Y Z]=chancenter(cell2mat({chans.X})', cell2mat({chans.Y})', cell2mat({chans.Z})',[]);
			 else
				 [X Y Z]=chancenter(cell2mat({chans.X})', cell2mat({chans.Y})', cell2mat({chans.Z})', extraargs{:});
			 end;
			 if isempty(X), return; end;
			 for index = 1:length(chans)
				 chans(index).X  = X(index);
				 chans(index).Y  = Y(index);
				 chans(index).Z  = Z(index);
			 end;
			case 'cart2topo',
			 [th rd]=cart2topo([cell2mat({chans.X})' cell2mat({chans.Y})' cell2mat({chans.Z})']);
			 if isempty(th), return; end;
			 for index = 1:length(chans)
				 chans(index).theta  = th(index);
				 chans(index).radius = rd(index);
			 end;
			case 'sph2topo',
			 disp('Warning: all radii considered to be one for this transformation');
			 try, [chan_num,angle,radius] = sph2topo([ones(length(chans),1)  cell2mat({chans.sph_phi})' cell2mat({chans.sph_theta})'], 1, 2); % using method 2
		     catch, error('Can not process empty values'); end;
			 for index = 1:length(chans)
				 chans(index).theta  = angle(index);
				 chans(index).radius = radius(index);
			 end;
			case 'topo2sph',
			 [sph_phi sph_theta] = topo2sph( [cell2mat({chans.theta})' cell2mat({chans.radius})'] );
			 for index = 1:length(chans)
				 chans(index).sph_theta  = sph_theta(index);
				 chans(index).sph_phi    = sph_phi  (index);
				 chans(index).sph_radius = 1;
			 end;
			case 'sph2cart',
			 [x y z] = sph2cart(cell2mat({chans.sph_theta})'/180*pi, cell2mat({chans.sph_phi})'/180*pi, cell2mat({chans.sph_radius})');
			 for index = 1:length(chans)
				 chans(index).X = x(index);
				 chans(index).Y = y(index);
				 chans(index).Z = z(index);
			 end;
			case 'cart2sph',
			 [th phi radius] = cart2sph(cell2mat({chans.X}), cell2mat({chans.Y}), cell2mat({chans.Z}));
			 for index = 1:length(chans)
				 chans(index).sph_theta     = th(index)/pi*180;
				 chans(index).sph_phi       = phi(index)/pi*180;
				 chans(index).sph_radius    = radius(index);
			 end;
		   end;
		  case 'transform'
		   tmpoper = args{curfield+1};
		   if isempty(tmpoper), return; end;
		   if iscell(tmpoper), tmpoper = tmpoper{1}; end;
		   if isempty(findstr(tmpoper, 'chans'))
			   try, X = cell2mat({chans.X}); catch, X(1:length(chans)) = NaN; end;
			   try, Y = cell2mat({chans.Y}); catch, Y(1:length(chans)) = NaN; end;
			   try, Z = cell2mat({chans.Z}); catch, Z(1:length(chans)) = NaN; end;
			   try, theta  = cell2mat({chans.theta}); catch, theta(1:length(chans)) = NaN; end;
			   try, radius = cell2mat({chans.radius}); catch, radius(1:length(chans)) = NaN; end;
			   try, sph_theta  = cell2mat({chans.sph_theta}); catch, sph_theta(1:length(chans)) = NaN; end;
			   try, sph_phi    = cell2mat({chans.sph_phi}); catch, sph_phi(1:length(chans)) = NaN; end;
			   try, sph_radius = cell2mat({chans.sph_radius}); catch, sph_radius(1:length(chans)) = NaN; end;
			   eval(tmpoper);
			   chans = struct('labels', { chans.labels }, 'X', mat2cell(X), 'Y', mat2cell(Y), 'Z', mat2cell(Z), ...
							  'theta', mat2cell(theta), 'radius', mat2cell(radius), 'sph_theta', mat2cell(sph_theta), ...
							  'sph_phi', mat2cell(sph_phi), 'sph_radius', mat2cell(sph_radius));
		   else 
			   eval(tmpoper);
		   end;
		  case 'shrink'
		   chans(1).shrink = args{ curfield+1 };		   
		  case 'delete'
		   chans(args{ curfield+1 })=[];
		  case 'changefield'
		   tmpargs = args{ curfield+1 };
		   if length( tmpargs ) < 3
			   error('pop_chanedit: not enough arguments to change field value');
		   end;
		   eval([ 'chans(' int2str(tmpargs{1}) ').'  tmpargs{2} '=' reformat(tmpargs{3} ) ';' ]);
		  case 'add'
		   tmpargs = args{ curfield+1 };
		   allfields = fieldnames(chans);
		   if length( tmpargs ) < length(allfields)+1
			   error('pop_chanedit: not enough arguments to change all field values');
		   end;
		   num = tmpargs{1};
		   chans(end+1) = chans(end);
		   chans(num+1:end) = chans(num:end-1);
		   for index = 1:length( allfields )
			   eval([ 'chans(' int2str(num) ').' allfields{index} '=' reformat(tmpargs{index+1}) ';' ]);
		   end;
		  case 'changechan'
		   tmpargs = args{ curfield+1 };
		   num = tmpargs{1};
		   allfields = fieldnames(chans);
		   if length( tmpargs ) < length(allfields)+1
			   error('pop_chanedit: not enough arguments to change all field values');
		   end;
		   for index = 1:length( allfields )
			   eval([ 'chans(' int2str(num) ').' allfields{index} '=' reformat(tmpargs{index+1}) ';' ]);
		   end;
		  case 'load'
		   tmpargs = args{ curfield+1 };
		   if ~isempty(tmpargs), 
			   if isstr(tmpargs)
				   chans = readlocs(tmpargs); 
			   else
				   if ~isempty(tmpargs{1}), 
					   chans = readlocs(tmpargs{:});
				   end;
			   end;
		   end;
		  case 'save'
		   tmpargs = args{ curfield+1 };
		   if isempty(tmpargs), return; end;
		   fid = fopen(tmpargs, 'w');
		   if fid ==-1, error('Can not open file'); end;
		   allfields = fieldnames(chans);
		   fprintf(fid, 'Number\t');
		   for field = 1:length(allfields)
			   fprintf(fid, '%s\t', allfields{field});
		   end;
		   fprintf(fid, '\n');
		   for index=1:length(chans)
			   fprintf(fid, '%d\t',  index);
			   for field = 1:length(allfields)
				   tmpval = getfield(chans, {index}, allfields{field});
				   if isstr(tmpval)
					   fprintf(fid, '%s\t',  tmpval);
				   else
					   fprintf(fid, '%3.3g\t',  tmpval);
				   end;
			   end;
			   fprintf(fid, '\n');
		   end;
		   if isempty(tmpargs), chantmp = readlocs(tmpargs); end;
		 end;
	 end;
end;
return;

function str = commandwarning(str);

% format the output field
% -----------------------
function strval = reformat( val )
	if isstr(val), strval = [ '''' val '''' ];
	else           strval = num2str(val);
	end;

function txt = inserttxt( txt, tokins, tokfind);
	locfind = findstr(txt, tokfind);
	for index = length(locfind):-1:1
		txt = [txt(1:locfind(index)-1) tokins txt(locfind(index):end)];
	end;
function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;
function [chans, shrinkfact]= checkchans(chans);
	shrinkfact = 'off';
	if isfield(chans, 'shrink'), 
		shrinkfact = chans(1).shrink; 
		chans = rmfield(chans, 'shrink');
	end;
    try, allfields = fieldnames(chans); catch, allfields = []; end;
	newstruct{1}  = 'labels';
	newstruct{3}  = 'theta';
	newstruct{5}  = 'radius';
	newstruct{7}  = 'X';
	newstruct{9}  = 'Y';
	newstruct{11} = 'Z';
	newstruct{13} = 'sph_theta';
	newstruct{15} = 'sph_phi';
	newstruct{17} = 'sph_radius';
	newstruct{18} = [];
	for index = 1:length(allfields)
		switch allfields{index}
		 case 'labels', pos=1;
		 case 'theta',  pos=3;
		 case 'radius', pos=5;
		 case 'X',      pos=7;
		 case 'Y',      pos=9;
		 case 'Z',      pos=11;
		 case 'sph_theta',    pos=13;
		 case 'sph_phi',      pos=15;
		 case 'sph_radius',   pos=17;
		 otherwise, disp(['Unrecognized field''' allfields{index} ''' in channel structure, removing it']); pos = -1;
		end;
		if pos ~= -1
			newstruct{pos+1} = eval(['{ chans.' allfields{index} '}']);
		end;
	end;
	chans = struct(newstruct{:});

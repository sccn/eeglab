% pop_chanedit() - Edit channel locations.
%
% Usage: >> newchans = pop_chanedit( chans, 'key1', value1, ...
%                                           'key2', value2, ... );
% Input:
%   chans - channel EEGLAB structure
%
% Optional inputs:
%   'convert'     - { conversion_type args } conversion type is 'xyz->polar'
%                   'sph2topo', 'topo2sph', 'sph2cart', 'cart2sph'. Args are
%                   only relevant for 'cart2topo', 'sph2topo', 'topo2sph' and
%                   can be found in the function having the same name.
%   'operation'   - string command for manipulating arrays. 'chan' is a full 
%                   electrode info. Fields can be manipulated using 'labels', 'theta'
%                   'radius' (polar angle and radius), 'X', 'Y', 'Z' (cartesian 3D) or
%                   'sph_theta', 'sph_phi', 'sph_radius' for spherical horizontal angle, 
%                   azimut and radius. Ex: 'chans(3) = chans(14)'. 'X = -X' or complex 
%                   formula separated by ';': 'TMP = X; X = Y; Y = TMP'
%   'changechan'  - {num value1 value2 value3 ...} Change the values of
%                   all fields the channel.
%   'changefield' - {num field value} change field value for channel number num. 
%                   (Ex: {34 'theta' 320.4}).
%   'add'         - {num label theta radius X Y Z sph_theta sph_phi sph_radius } 
%                   Insert channel before channel number num with the specified values. 
%                   If the number of values if less than ten, fields are filled with 0.
%   'delete'      - vector of indices of channel to delete
%   'load'        - { textfile 'format' } format can be 'loc' for polar coordinate
%                   file, 'elp' for Polhemus text files, 'sph1' for spherical files
%                   with no radius.
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
	 otherwise, error(['Unrecognized field''' allfields{index} ''' in input structure']);
	end;
	newstruct{pos+1} = eval(['{ chans.' allfields{index} '}']);
end;
chans = struct(newstruct{:});

allfields = fieldnames(chans);
if nargin < 2
	totaluserdat = {};
	ingui = 1;
	guimodif = 0;
	while ingui
		commentfields = { 'Channel label', 'Polar angle (theta)', 'Polar radius (radius)', 'Cartesian X (X)', 'Cartesian Y (Y)', 'Cartesian Z (Z)', ...
						  'Spherical horiz. angle (sph_theta)', 'Spherical azimuth angle (sph_phi)', 'Spherical radius (sph_radius)' };
		% transfer channel to global workspace
		global chantmp;
		chantmp = chans;
		evalin('base', [ 'global chantmp ;' ]);
		
		% add field values
		% ----------------
		geometry = { 1 };
		tmpstr = sprintf('Edit channels info');
		uilist = { { 'Style', 'text', 'string', tmpstr, 'fontweight', 'bold'  } };
		endgui = 'set(findobj(''parent'', gcbf, ''tag'', ''ok''), ''userdata'', ''stop'');';
		operation = ['inputdlg({strvcat(''Enter operation (see help on previous page)''' ...
					 ',''(Ex: TMP=X; X=-Y; Y=TMP or Y(3) = X(2)'')}, ''Operation'', 1, { '''' });'];
		uiconvert = { { 'Style', 'pushbutton', 'string', 'xyz->polar', 'callback', ...
						['comtmp = {''convert'' {''cart2topo'' ''gui'' ''on''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'sph->polar', 'callback', ...
						['comtmp = {''convert'' {''sph2topo'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'polar->sph', 'callback', ...
						['cmptmp = {''convert'' {''topo2sph'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'sph->xyz'  , 'callback', ...
						['comtmp = {''convert'' {''sph2cart'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'xyz->sph'  , 'callback', ...
						['comtmp = {''convert'' {''cart2sph'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'operation ', 'callback', ...
						['comtmp = {''operation'' ' operation '};' endgui] }, ...
					  {} { } { } };
		%{ 'Style', 'pushbutton', 'string', 'UNDO LAST ', 'callback', '' } { } { } };
		uiconvert{1}{6}		  
		for index = 1:length(allfields)
			geometry = { geometry{:} [1.5 1 0.2 1] };
			uilist   = { uilist{:}, ...
						 { 'Style', 'text', 'string', commentfields{index} }, ...
						 { 'Style', 'edit', 'tag', allfields{index}, 'string', num2str(getfield(chans,{1}, allfields{index})), 'callback', ...
						   [ 'valnum   = str2num(get(findobj(''parent'', gcbf, ''tag'', ''numval''), ''string''));' ...
							 'editval = get(gcbo, ''string'');' ...
							 'if ~isempty(str2num(editval)), editval =str2num(editval);  end;' ...
							 'eval([ ''chantmp(valnum).' allfields{index} '= editval;'']);' ...
							 'olduserdat = get( gcbf, ''userdata'');' ...
							 'if isempty(olduserdat), olduserdat = {}; end;' ...
							 'set( gcbf, ''userdata'', { olduserdat{:} ''changefield'' { valnum ''' allfields{index} ''' editval }});' ...
							 'clear editval valnum olduserdat;' ] } { } uiconvert{index} };
		end;
		
		% add buttons
		% -----------
		geometry =  { geometry{:} [1] [1 0.7 0.4 1.8 0.2 0.7 1] [1 0.7 0.7 1 0.7 0.7 1] };
		callpart1 = [ 'valnum   = str2num(get(findobj(''parent'', gcf, ''tag'', ''numval''), ''string''));' ];
		callpart2 = [ 'set(findobj(''parent'', gcf, ''tag'', ''numval''), ''string'', num2str(valnum));' ];
		for index = 1:length(allfields)
			callpart2 = [ callpart2  'set(findobj(''parent'', gcf, ''tag'', ''' allfields{index} ...
						  '''), ''string'', num2str(chantmp(valnum).' allfields{index} '));' ];
		end;
		callpart2 = [ callpart2 'set(findobj(''parent'', gcf, ''tag'', ''scantitle''), ' ...
					  '''string'', [''Channel Num (out of '' int2str(length(chantmp)) '')'']);' ]; 
		callpart2 = [ callpart2 'clear valnum;' ];
		
		uilist   = { uilist{:}, ...
					 { }, ...
					 { },{ },{ }, {'Style', 'text', 'string', ['Channel Num (out of ' int2str(length(chans)) ')'], 'fontweight', 'bold', 'tag', 'scantitle' }, { },{ },{ }, ...
					 { 'Style', 'pushbutton', 'string', 'Delete chan',  'callback', [callpart1 'chantmp(valnum) = []; valnum = min(valnum,length(chantmp));' ...
					'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
					'set( gcbf, ''userdata'', { olduserdat{:} ''delete'', valnum }); clear olduserdat' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '<<', 'callback', [callpart1 'valnum = max(valnum-10,1);' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '<',  'callback', [callpart1 'valnum = max(valnum-1,1);' callpart2 ] }, ...
					 { 'Style', 'edit', 'string', '1', 'tag', 'numval', 'callback', [callpart1 'valnum = min(str2num(get(gcbo, ''string'')),length(chantmp));' callpart2 ] }, ...
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
		plot2dcom = [ 'tmpshrink = ''off'';' ...
					  'if get(findobj(''parent'', gcbf,  ''string'', ''Auto shrink''), ''value'')' ...
					  '   tmpshrink = ''force'';' ...
					  'else ' ...
					  '   if ~isempty(get(findobj(''parent'', gcbf, ''tag'', ''shrinkfactor''), ''string''))' ...
					  '       tmpshrink = str2num(get(findobj(''parent'', gcbf, ''tag'', ''shrinkfactor''), ''string''));' ...
					  '   end;' ...
					  'end;' ...
					  'figure; topoplot([],chantmp, ''style'', ''blank'', ''electrodes'', ''labelpoint'', ''shrink'', tmpshrink); tmpshrink, clear tmpshrink;'];
		geometry = { geometry{:} [1] [1 1.4 0.6 1 1] [1] [1 1 1 1]};
		uilist = {  uilist{:},...
					{ } ...
					{ 'Style', 'pushbutton', 'string', 'Plot 2D', 'callback', plot2dcom },... 
					{ 'Style', 'text', 'string', '2D shink factor (-1 to 1)'} ...
					{'Style', 'edit', 'string', '', 'tag', 'shrinkfactor'} ...
					{'Style', 'checkbox', 'string', 'Auto shrink'} ...
					{ 'Style', 'pushbutton', 'string', 'Plot 3D', 'callback', 'plotchans3d([cell2mat({chantmp.X})'' cell2mat({chantmp.Y})'' cell2mat({chantmp.Z})''], {chantmp.labels});' } ...
					{}, { 'Style', 'pushbutton', 'string', 'Load', 'callback', ...
						  ['[tmpf tmpp] = uigetfile(''*'', ''Load a channel location file'');' ...
						   'comtmp = {''load'' [tmpp tmpf] };' endgui ] }, ...
					{ 'Style', 'pushbutton', 'string', 'Load help', 'callback', 'pophelp(''readlocs.m'');' }, ...	
					{ 'Style', 'pushbutton', 'string', 'Save ASCII', 'callback', ...	
					  ['[tmpf tmpp] = uiputfile(''*'', ''Save channel location file'');' ...
					   'comtmp = {''save'' [tmpp tmpf] };' endgui ] }, ...
					{ 'Style', 'pushbutton', 'enable', 'off', 'string', 'Save BESA' }, ...	
				 };
		
		if guimodif
			set(gcf, 'userdat', {});
			eval([ callpart1 callpart2 ]);
			[results userdat returnmode] = inputgui( geometry, uilist, 'pophelp(''pop_chanedit'');', 'Edit channel info -- pop_chanedit()', {}, gcf );
		else 
			[results userdat returnmode] = inputgui( geometry, uilist, 'pophelp(''pop_chanedit'');', 'Edit channel info -- pop_chanedit()', {}, 'noclose' );
		end; 
		if length(results) == 0, return; end;
    
		% transfer events back from global workspace
		chans = evalin('base', 'chantmp;');
		%evalin('base', 'clear chantmp;');
		totaluserdat = { totaluserdat{:} userdat{:}};
		if evalin('base', 'exist(''comtmp'')') == 1
			tmpcom = evalin('base', 'comtmp;');
			chans = pop_chanedit(chans, tmpcom{:}); % apply modification to channel structure
			totaluserdat = { totaluserdat{:} tmpcom };
			evalin('base', 'clear comtmp;');
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
	if ~isempty(findobj('parent', gcf, 'tag','shrinkfactor'))
		close(gcf);
	end;
	if ~isempty( totaluserdat )
		com = sprintf('%s=pop_chanedit(%s, %s);', vararg2str(totaluserdat));
	end;
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
		   if isstr(extraargs{1}) & strcmp(extraargs{1}, 'gui') & ~strcmp(method, 'cart2topo')
			   tmpButtonName=questdlg( ['This will modify fields in the channel structure' 10 ...
					'Are you sure, you want to apply that function ?'], 'Confirmation', 'Cancel', 'Yes','Yes');
			   if ~strcmp(tmpButtonName, 'Yes'), return; end;
		   end;
		   switch method
			case 'cart2topo',
			 if isempty(extraargs)
				 [th rd]=cart2topo([cell2mat({chans.X})' cell2mat({chans.Y})' cell2mat({chans.Z})']);
			 else
				 [th rd]=cart2topo([cell2mat({chans.X})' cell2mat({chans.Y})' cell2mat({chans.Z})'], extraargs{:});
			 end;
			 if isempty(th), return; end;
			 for index = 1:length(chans)
				 chans(index).theta  = th(index);
				 chans(index).radius = rd(index);
			 end;
			case 'sph2topo',
			 disp('Warning: all radii considered to be one for this transformation');
			 try, [chan_num,angle,radius] = sph2topo([ones(length(chans),1) cell2mat({chans.sph_theta})' cell2mat({chans.sph_phi})'], 1, 2); % using method 2
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
		  case 'operation'
		   tmpoper = args{curfield+1};
		   if iscell(tmpoper), tmpoper = tmpoper{1}; end;
		   if isempty(findstr(tmpoper, 'chans'))
			   try, X = cell2mat({chans.X}); catch, X(1:length(chans)) = NaN; end;
			   try, Y = cell2mat({chans.Y}); catch, Y(1:length(chans)) = NaN; end;
			   try, Z = cell2mat({chans.Z}); catch, Z(1:length(chans)) = NaN; end;
			   try, theta  = cell2mat({chans.theta}); catch, theta(1:length(chans)) = NaN; end;
			   try, radius = cell2mat({chans.radius}); catch, radius(1:length(chans)) = NaN; end;
			   try, sph_theta  = cell2mat({chans.sph_theta}); catch, sph_theta(1:length(chans)) = NaN; end;
			   try, sph_phi    = cell2mat({chans.sph_theta}); catch, sph_phi(1:length(chans)) = NaN; end;
			   try, sph_radius = cell2mat({chans.sph_theta}); catch, sph_radius(1:length(chans)) = NaN; end;
			   eval(tmpoper);
			   chans = struct('labels', { chans.labels }, 'X', mat2cell(X), 'Y', mat2cell(Y), 'Z', mat2cell(Z), ...
							  'theta', mat2cell(theta), 'radius', mat2cell(radius), 'sph_theta', mat2cell(sph_theta), ...
							  'sph_phi', mat2cell(sph_phi), 'sph_radius', mat2cell(sph_radius));
		   else 
			   eval(tmpoper);
		   end;
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
		   if ~isempty(tmpargs), chans = readlocs(tmpargs); end;
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
	 ButtonName=questdlg( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;
	
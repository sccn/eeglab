% pop_chanedit() - Edit channel locations (chanlocs) structure of an EEGLAB dataset.
%                  For channel location structure and file format, 
%                  see >> help readlocs  % help message of readlocs() 
%
% Usage: >> newchans = pop_chanedit( EEG, 'key1', value1, ...
%                        'key2', value2, ... ); % edit dataset containing chanlocs
%        >> newchans = pop_chanedit( chanlocs, 'key1', value1, ...
%                        'key2', value2, ... ); % edit separate chanlocs struct
% Graphic interface:
%   "Channel information ('field name')" - [edit boxes] display channel field 
%                   content of the current channel. Use 'transform' from the command
%                   line to modify these fields.
%   "Opt. 3D center" - [button] recenter 3-D channel coordinates. Uses the chancenter()
%                 function. Command line equivalent is 'convert', { 'chancenter' [xc
%                 yc zc] }, [xc yc zc] being the center of the sphere (use empty to
%                 find the center so that electrodes best match a sphere).
%   "rotate axis" - [button] force one electrode to one position and rotate other
%                 electrodes accordingly. Command line equivalent is 'forcelocs'.
%   "Transform axis" - [button] perform any operation on channel fields. Command
%                 line equivalent is 'transform'.
%   "xyz->polar & sph." - [button] convert 3-D cartesian coordinates to polar and
%                 3-D spherical coordinates. This is usefull when you edit the 
%                 coordinates manually. Command line equivalent is 'convert', 
%                 'cart2all'.
%   "sph.->polar & xyz" - [button] convert 3-D spherical coordinates to polar and
%                 3-D cartesian coordinates. Command line equivalent is 'convert', 
%                 'sph2all'.
%   "polar->sph & xyz" - [button] convert 2-D polar coordinates to 3-D spherical and
%                 3-D cartesian coordinates. Command line equivalent is 'convert', 
%                 'topo2all'. Note that if spherical radii are absent, they are forced
%                 to 1.
%   "Delete chan" - [button] delete channel. Command line equivalent: 'delete'.
%   "Insert chan" - [button] insert channel before current channel. Command line 
%                 equivalent: 'insert'.
%   "<<" - [button] scroll channel backward (-10).
%   "<" - [button] scroll channel backward (-1).
%   ">" - [button] scroll channel forward (+1).
%   ">>" - [button] scroll channel forward (+10).
%   "Append chan" - [button] append channel after the current channel. Command line 
%                 equivalent: 'append'.
%   "Plot 2D" - [button] plot channel positions in 2-D using topoplot() function.  
%   "2-D shrink factor" - [edit box] enter a value to shrink channel location in the
%                 2-D polar view. This DOES NOT AFFECT channel locations and is only
%                 used for visualization. This parameter is attached to the dataset
%                 and is then used in all 2-D scalp plots. Command line equivalent
%                 is 'shink'.
%   "Auto shrink" - [button] automatically computes the shrinking factor so all
%                 channels are visible on the 2-D topographical plot.
%   "Plot 3D" - [button] plot channel positions in 3-D using plotchans3d() function.                 
%   "Read locations" - [button] read location file using readlocs() function. 
%                 Command line equivalent is 'load'.
%   "Read help" - [button] readlocs() function help.
%   "Save .ced" - [button] save channel location as ".ced" format which is the native
%                 EEGLAB format. Command line equivalent is 'save'.
%   "Save others" - [button] save channel location as other formats using the 
%                 pop_writelocs() function (see also readlocs() for channel formats).
%   "Cancel" - [button] cancel all editing.
%   "Help" - [button] this help message.
%   "OK" - [button] save editing and propagate to parent.
% 
% Input:
%   EEG      - EEG dataset
%   chanlocs - EEG.chanlocs structure
%
% Optional inputs:
%   'convert'     - {conversion_type [args]} Conversion type may be: 'cart2topo'
%                   'sph2topo', 'topo2sph', 'sph2cart', 'cart2sph', or 'chancenter'. 
%                   See help messages for these functions. Args are only relevant 
%                   for 'chancenter'. See also graphical interface button for more
%                   info.
%   'transform'   - String command for manipulating arrays. 'chan' is full channel 
%                   info. Fields that can be manipulated are 'labels', 'theta'
%                   'radius' (polar angle and radius), 'X', 'Y', 'Z' (cartesian 
%                   3-D) or 'sph_theta', 'sph_phi', 'sph_radius' for spherical 
%                   horizontal angle, azimuth and radius. 
%                   Ex: 'chans(3) = chans(14)', 'X = -X' or a multi-step transform
%                   with steps separated by ';': Ex. 'TMP = X; X = Y; Y = TMP'
%   'changechan'  - {num value1 value2 value3 ...} Change the values of all fields 
%                   for the given channel num: mimimally {num label theta radius}.
%                   Ex: 'changechan' {12 'PXz' -90 0.30}
%   'changefield' - {num field value} Change field value for channel number num. 
%                   Ex: {34 'theta' 320.4}.
%   'insert'      - {num label theta radius X Y Z sph_theta sph_phi sph_radius } 
%                   Insert new channel before channel number num with the specified 
%                   values. If the number of values if less than 10, remaining 
%                   fields are 0. (previously this parameter was termed 'add').
%   'append'      - {num label theta radius X Y Z sph_theta sph_phi sph_radius } 
%                   same as 'insert' but add the the new channel after the current
%                   channel.
%   'delete'      - [int_vector] Vector of channel numbers to delete.
%   'forcelocs'   - [cell] call forcelocs() function to force a particular channel
%                   to be at a particular location on the sphere (and rotate other
%                   channels accordingly).
%   'shrink'      - Topographical polar shrink factor (see >> help topoplot) 
%   'load'        - [filename|{filename, 'key', 'val'}] Load channel location file
%                   optional arguments (such as file format) to the function 
%                   readlocs() can be specified if the input is a cell array.
%   'save'        - 'filename' Save text file with channel info.
%   'eval'        - [string] evaluate string ('chantmp' is the name of the channel
%                   location structure.
%   'lookup'      - [integer array] lookup channel indices standard location. []
%                   looks up all channels.
%
% Outputs:
%   newchans      - new EEGLAB channel locations structure
%
% Ex:  EEG = pop_chanedit(EEG,'load', { 'dummy.elp' 'elp' }, 'delete', [3 4], ...
%          'convert', { 'xyz->polar' [] -1 1 }, 'save', 'mychans.loc' )
%        % Load polhemus file, delete two channels, convert to polar (see 
%        % cart2topo() for arguments) and save into 'mychans.loc'.
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
% Revision 1.88  2004/02/13 02:36:21  arno
% updating location after chancenter
%
% Revision 1.87  2004/02/13 00:10:31  arno
% implementing new history for channel recentering
%
% Revision 1.86  2004/02/12 23:57:47  arno
% *** empty log message ***
%
% Revision 1.85  2004/01/29 00:26:14  arno
% text in message
%
% Revision 1.84  2004/01/01 19:18:50  scott
% same
%
% Revision 1.83  2004/01/01 19:16:09  scott
% 3d center -> Opt. 3D center
%
% Revision 1.82  2003/12/20 01:47:09  arno
% debug call with empty array
%
% Revision 1.81  2003/12/18 01:21:24  arno
% same
%
% Revision 1.80  2003/12/18 01:17:23  arno
% warning message
%
% Revision 1.79  2003/12/18 01:10:37  arno
% debug text
%
% Revision 1.78  2003/12/18 01:08:20  arno
% channel location lookup
%
% Revision 1.77  2003/12/17 01:35:58  arno
% debug plot3d for empty chans
%
% Revision 1.76  2003/12/12 01:17:22  arno
% nothing
%
% Revision 1.75  2003/12/10 03:19:56  arno
% graphical interface help
%
% Revision 1.74  2003/12/05 23:41:34  arno
% debug convert2all
%
% Revision 1.73  2003/12/05 23:39:09  arno
% coordinate conversion debug
%
% Revision 1.72  2003/12/05 22:58:23  arno
% same thing
%
% Revision 1.71  2003/12/05 22:54:08  arno
% debug transform
%
% Revision 1.70  2003/12/05 18:25:45  arno
% ?
%
% Revision 1.69  2003/12/05 18:20:21  arno
% same thing
%
% Revision 1.68  2003/12/05 18:18:52  arno
% checkchans problem
%
% Revision 1.67  2003/12/02 22:43:38  arno
% forcing compatibility with urchan
%
% Revision 1.66  2003/12/02 19:25:38  arno
% debug readlocs and chanedit
%
% Revision 1.65  2003/12/02 18:02:06  arno
% chancenter history
%
% Revision 1.64  2003/12/02 17:55:37  arno
% verbose on for convertlocs
%
% Revision 1.63  2003/12/02 17:54:34  arno
% conversion of coordinates
%
% Revision 1.62  2003/12/02 17:15:08  arno
% adding insert & append button
%
% Revision 1.61  2003/12/02 03:31:40  arno
% current figure problem
%
% Revision 1.60  2003/12/02 03:21:26  arno
% better gui for readlocs
%
% Revision 1.59  2003/10/16 23:39:45  arno
% eeglabsources
% OA
% eeglabsources
%
% Revision 1.58  2003/08/08 17:00:37  arno
% fix header
%
% Revision 1.57  2003/08/05 18:30:46  arno
% nothing
%
% Revision 1.56  2003/08/04 18:49:56  arno
% automatic conversion for 'transform' parameter
%
% Revision 1.55  2003/07/19 01:24:47  arno
% remving double &
%
% Revision 1.54  2003/07/16 18:45:50  arno
% debug string
% []
%
% Revision 1.53  2003/07/16 16:45:00  arno
% debug readlocs pop up window
%
% Revision 1.52  2003/07/16 01:38:07  scott
% help topoplot
% Revision 1.51  2003/07/10 18:30:03  arno
% debuging rotate axis cancelation
%
% Revision 1.50  2003/05/13 23:49:13  arno
% [Aallowing to write channel location file
%
% Revision 1.49  2003/05/13 16:55:35  arno
% debug shrink factor
%
% Revision 1.48  2003/04/17 01:55:35  arno
% automatically convert coordinates for 3-d center
%
% Revision 1.47  2003/04/16 02:06:08  arno
% adding forcelocs option
%
% Revision 1.46  2003/04/10 17:29:47  arno
% header edit
%
% Revision 1.45  2003/03/06 00:36:34  arno
% further checking
%
% Revision 1.44  2003/03/06 00:35:30  arno
% same
%
% Revision 1.43  2003/03/06 00:34:28  arno
% handling numeric shrink factors
%
% Revision 1.42  2003/02/23 08:23:08  scott
% header edit -sm
%
% Revision 1.41  2003/01/03 22:43:35  arno
% removing error message
% ,
%
% Revision 1.40  2003/01/02 17:31:49  arno
% text in pop-up window for reading channel location
%
% Revision 1.39  2002/12/27 22:58:07  arno
% removing debugging message
%
% Revision 1.38  2002/11/15 02:32:27  arno
% debug for command line call
%
% Revision 1.37  2002/11/15 01:38:14  scott
% Can not -> cannot
%
% Revision 1.36  2002/11/13 17:44:11  arno
% header edition -sm, ad
%
% Revision 1.35  2002/11/13 17:12:50  scott
% help msg - changechan
%
% Revision 1.34  2002/11/13 16:57:24  scott
% help msg
%
% Revision 1.33  2002/11/13 14:57:54  scott
% help msg edit
%
% Revision 1.32  2002/11/12 23:02:21  arno
% debugging message when number of channel does not match
%
% Revision 1.31  2002/10/23 15:51:35  arno
% more bug fixed
%
% Revision 1.30  2002/10/23 15:21:23  arno
% saving error fixed
%
% Revision 1.29  2002/10/16 16:17:13  arno
% debug command line call
%
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
allfields = { 'labels' 'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' };
[chans shrinkfact]= checkchans(chans, allfields);

if nargin < 2
    
	totaluserdat = {};
    % lookup channel locations if necessary
    % -------------------------------------
    if ~all(cellfun('isempty', {chans.labels})) & all(cellfun('isempty', {chans.theta}))
        standardchans = { 'Fp1' 'Fpz' 'Fp2' 'Nz' 'AF9' 'AF7' 'AF3' 'AFz' 'AF4' 'AF8' 'AF10' 'F9' 'F7' 'F5' ...
                          'F3' 'F1' 'Fz' 'F2' 'F4' 'F6' 'F8' 'F10' 'FT9' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' ...
                          'FC4' 'FC6' 'FT8' 'FT10' 'T9' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'T10' ...
                          'TP9' 'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10' 'P9' 'P7' 'P5' ...
                          'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO9' 'PO7' 'PO3' 'POz' 'PO4' 'PO8' 'PO10' ...
                          'O1' 'Oz' 'O2' 'O9' 'O10' 'CB1' 'CB2' 'Iz' };
        [tmp1 ind1 ind2] = intersect( lower(standardchans), lower({ chans.labels }));
        if ~isempty(tmp1)
            res = questdlg2(strvcat('Only channel labels are currenlty present but some', ...
                                   'of these labels have standart locations. Do you want to look up', ...
                                   'standard coordinates for these channels?'), 'Look up channel locations?', ...
                           'No', 'Yes', 'Yes');
            if strcmpi(res, 'yes')
                chans = pop_chanedit(chans, 'lookup', []);	
                totaluserdat = { 'lookup' [] };
            end;
        end;
    end;
    
	ingui = 1;
	guimodif = 0;
	while ingui
		commentfields = { 'Channel label ("label")', 'Polar angle ("theta")', 'Polar radius ("radius")', ...
                          'Cartesian X ("X")', 'Cartesian Y ("Y")', 'Cartesian Z ("Z")', ...
						  'Spherical horiz. angle ("sph_theta")', 'Spherical azimuth angle ("sph_phi")', ...
                          'Spherical radius ("sph_radius")' };
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
		transform = [ 'inputdlg2({''Enter transform: (Ex: TMP=X; X=-Y; Y=TMP or Y(3) = X(2), etc.'' }, ' ...
                      '''Transform'', 1, { '''' }, ''pop_chanedit'');'];
        guicenter = [ '[chantmp newcenter tmpcom] = pop_chancenter(chantmp);' ...
                      'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
                      'set( gcbf, ''userdata'', { olduserdat{:} ''eval'', tmpcom });' ...
                      'clear olduserdat tmpcom newcenter; comtmp = {};' endgui ];
        %guicenter = [ '[tmpX tmpY tmpZ newcenter tmpoptim]=chancenter(cell2mat({chanstmp.X})'', ' ...
        %              '               cell2mat({chanstmp.Y})'', cell2mat({chanstmp.Z})'',[]);' ...
        %              'chanstmp = pop_chanedit(chanstmp, ''convert'', { ''chancenter'' newcenter 0 });' ...
        %              'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
        %              'if tmpoptim, newcenter = []; end;' ...
        %              'set( gcbf, ''userdata'', { olduserdat{:} ''convert'', { ''chancenter'' newcenter 0 } });' ...
        %              'clear tmpcell olduserdat tmpoptim tpmpX tmpY tmpZ newcenter;' ];

		uiconvert = { { 'Style', 'pushbutton', 'string', 'Opt. 3D center', 'callback', ...
						guicenter }, ...
					  { 'Style', 'pushbutton', 'string', 'Rotate axis'  , 'callback', ...
						['[ comtmp tmpforce ] = forcelocs(chantmp); if ~isempty(tmpforce), ' ...
                         'comtmp = {''forcelocs'' tmpforce{1} }; else comtmp = {''forcelocs'' [] }; end; clear tmpforce;' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'Transform axes', 'callback', ...
						['comtmp = {''transform'' ' transform '};' endgui] }, ...
					  {  }, ...
					  { 'Style', 'pushbutton', 'string', 'xyz -> polar & sph.', 'callback', ...
						['comtmp = {''convert'' {''cart2all'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'sph. -> polar & xyz', 'callback', ...
						['comtmp = {''convert'' {''sph2all'' ''gui''}};' endgui] }, ...
					  { 'Style', 'pushbutton', 'string', 'polar -> sph. & xyz', 'callback', ...
						['comtmp = {''convert'' {''topo2all'' ''gui''}};' endgui] }, ...
					  {  }, ...
					  { } { } };
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
		geometry =  { geometry{:} [1] [1.15 0.5 0.6 1.9 0.4 0.4 1.15] [1.15 0.7 0.7 1 0.7 0.7 1.15] };
		callpart1 = [ 'valnum   = str2num(char(get(findobj(''tag'', ''chaneditnumval''), ''string'')));' ];
		callpart2 = [ 'set(findobj(''tag'', ''chaneditnumval''), ''string'', num2str(valnum));' ];
		for index = 1:length(allfields)
			callpart2 = [ callpart2  'set(findobj(''tag'', ''chanedit' allfields{index} ...
						  '''), ''string'', num2str(chantmp(valnum).' allfields{index} '));' ];
		end;
		callpart2 = [ callpart2 'set(findobj(''tag'', ''chaneditscantitle''), ' ...
					  '''string'', [''Channel number (of '' int2str(length(chantmp)) '')'']);' ...
                      'set(findobj(''tag'', ''chaneditnumval''), ''string'', int2str(valnum));'  ]; 
		callpart2 = [ callpart2 'clear orivalnum valnum;' ];
		cb_del    = [ callpart1 'orivalnum = valnum; chantmp(valnum) = [];' ...
					'valnum = min(valnum,length(chantmp));' ...
                    'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
					'set( gcbf, ''userdata'', { olduserdat{:} ''delete'', orivalnum }); clear olduserdat' callpart2 ];
        cb_insert = [callpart1 ...
					'chantmp(end+3) = chantmp(end);' ...
					'chantmp(valnum+1:end-2) = chantmp(valnum:end-3);' ...
					'chantmp(valnum) = chantmp(end-1);' ...
					'chantmp = chantmp(1:end-2);' ...
					'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
					'tmpcell = cell(1,1+length(fieldnames(chantmp))); tmpcell{1} =valnum;' ...
					'set( gcbf, ''userdata'', { olduserdat{:} ''insert'', tmpcell }); clear tmpcell olduserdat' callpart2 ];
        cb_append = [callpart1 ...
					'chantmp(end+3) = chantmp(end);' ...
					'chantmp(valnum+2:end-2) = chantmp(valnum+1:end-3);' ...
					'chantmp(valnum+1) = chantmp(end-1);' ...
					'chantmp = chantmp(1:end-2);' ...
					'olduserdat = get( gcbf, ''userdata''); if isempty(olduserdat), olduserdat = {}; end;' ...
					'tmpcell = cell(1,1+length(fieldnames(chantmp))); tmpcell{1} =valnum;' ...
					'set( gcbf, ''userdata'', { olduserdat{:} ''append'', tmpcell }); valnum = valnum+1; clear tmpcell olduserdat' callpart2 ];
        
		uilist   = { uilist{:}, ...
					 { }, ...
					 { 'Style', 'pushbutton', 'string', 'Delete chan',  'callback', cb_del }, ...
					 { },{ }, ...
                     { 'Style', 'text'      , 'string', ['Channel number (of ' int2str(length(chans)) ')'], ...
					                                     'fontweight', 'bold', 'tag', 'chaneditscantitle' }, { },{ },{ }, ...
					 { 'Style', 'pushbutton', 'string', 'Insert chan',  'callback', cb_insert } ...
					 { 'Style', 'pushbutton', 'string', '<<', 'callback', [callpart1 'valnum = max(valnum-10,1);' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '<',  'callback', [callpart1 'valnum = max(valnum-1,1);' callpart2 ] }, ...
					 { 'Style', 'edit'      , 'string', '1', 'tag', 'chaneditnumval', 'callback', ...
					                                [callpart1 'valnum = min(str2num(get(gcbo, ''string'')),length(chantmp));' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '>',  'callback', [callpart1 'valnum = min(valnum+1,length(chantmp));' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', '>>', 'callback', [callpart1 'valnum = min(valnum+10,length(chantmp));' callpart2 ] }, ...
					 { 'Style', 'pushbutton', 'string', 'Append chan',  'callback', cb_append }, ...
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
        guireadlocs = [ '[tmpf tmpp] = uigetfile(''*'', ''Load a channel location file''); drawnow;' ...
                        'if ~isequal(tmpf, 0),' ...
                        '   tmpformats = { ''autodetect'', ''loc eeglab'', ''sph eeglab'', ''sfp egi/besa'', ''xyz eeglab'', ''asc neuroscan'', ' ...
                        '         ''dat neuroscan'', ''elc (eetrak)'', ''polhemusx (elp)'', ''polhemusy (elp)'', ''besa (elp)'', ''chanedit (ced)'' };' ...
                        '   tmpfmt = inputgui({[1 1]}, { { ''style'', ''text'', ''string'', [ ''File format'' 10 10 10 10 ] },' ...
                        '            { ''style'', ''listbox'', ''string'', strvcat(tmpformats) }}, ''pophelp(''''readlocs'''')'',' ...
                        '            ''Read electrode file'', [], ''normal'', 4);' ...
						'   if isempty(tmpfmt), comtmp = {''load'' [] };' ...
						'   else comtmp = {''load'' { [tmpp tmpf ] ''filetype'' tmpformats{tmpfmt{1}} } };' ... 
                        '   end; clear tmpfmt tmpp tmpf tmpformats;' ...
                        'else comtmp = {''load'' [] };' ...
                        'end;' endgui ];
        plot3d =     [ 'tmpind = find(~cellfun(''isempty'', { chantmp.X }));' ...
					   'if ~isempty(tmpind),' ...
                       '   plotchans3d([cell2mat({chantmp(tmpind).X})'' cell2mat({chantmp(tmpind).Y})'' cell2mat({chantmp(tmpind).Z})''],' ...
					                   '{chantmp(tmpind).labels}); else disp(''cannot plot: no XYZ coordinates'');' ...
                       'end;' ];
		uilist = {  uilist{:},...
					{ } ...
					{ 'Style', 'pushbutton', 'string', 'Plot 2D', 'callback', plot2dcom },... 
					{ 'Style', 'text', 'string', '2D shink factor (-1 to 1)'} ...
					{ 'Style', 'edit', 'string', fastif(isnumeric(shrinkfact), num2str(shrinkfact), ''), ...
					  'tag', 'shrinkfactor', 'callback', 'tmpshrink=get(gcbo, ''string''); set(findobj(''tag'', ''shrinkbut''), ''value'', 0);' } ...
					{ 'Style', 'checkbox', 'tag', 'shrinkbut', 'string', 'Auto shrink', 'value', strcmp(shrinkfact, 'force'), ...
					  'callback', 'if get(gcbo, ''value''), tmpshrink=''force''; else tmpshrink=''off''; end; set(findobj(''tag'', ''shrinkfactor''), ''string'', '''');' } ...
					{ 'Style', 'pushbutton', 'string', 'Plot 3D (XYZ)', 'callback', plot3d } ...
					{}, { 'Style', 'pushbutton', 'string', 'Read locations', 'callback', guireadlocs }, ...
					{ 'Style', 'pushbutton', 'string', 'Read help', 'callback', 'pophelp(''readlocs.m'');' }, ...	
					{ 'Style', 'pushbutton', 'string', 'Save .ced', 'callback', ...	
					  ['[tmpf tmpp] = uiputfile(''*.ced'', ''Save channel locs in EEGLAB .ced format''); drawnow;' ...
					   'comtmp = {''save'' [tmpp tmpf] }; clear tmpfmt tmpp tmpf;' endgui ] }, ...
					{ 'Style', 'pushbutton', 'string', 'Save others' 'callback', ...	
					  ['com = pop_writelocs(chantmp); comtmp = {''eval'' com }; clear tmpfmt tmpp tmpf;' endgui ] }, ...	
				 };
		
%						  ['[tmpf tmpp] = uigetfile(''*'', ''Load a channel location file''); drawnow;' ...
%						   'tmpfmt = inputdlg2({[''File format [loc|sph|sfp|xyz|asc|dat|polhemusx|polhemusy|besa|chanedit] ([]=use extension)'' ]}, ''File format'', 1, { '''' }, ''readlocs'');' ...
%						   'if isempty(tmpfmt), tmpfmt = {[]}; tmpp = []; tmpf = []; end;' ...
%						   'comtmp = {''load'' { [tmpp tmpf ] ''filetype'' tmpfmt{1} } }; clear tmpfmt tmpp tmpf;' endgui ] }, ...
		if guimodif
			set(currentfig, 'userdat', {});
			eval([ callpart1 callpart2 ]);
			[results userdat returnmode] = inputgui( geometry, uilist, 'pophelp(''pop_chanedit'');', 'Edit channel info -- pop_chanedit()', {}, currentfig );
		else 
			[results userdat returnmode] = inputgui( geometry, uilist, 'pophelp(''pop_chanedit'');', 'Edit channel info -- pop_chanedit()', {}, 'noclose' );
			if ~isempty(get(0, 'currentfigure')) currentfig = gcf; end;
		end; 
		if length(results) == 0, return; end;
		
		% transfer events back from global workspace
		chans = evalin('base', 'chantmp');
		%evalin('base', 'clear chantmp');
		totaluserdat = { totaluserdat{:} userdat{:}};
		if evalin('base', 'exist(''comtmp'')') == 1
			tmpcom = evalin('base', 'comtmp');
            if ~isempty(tmpcom)
                try, 
                    chans = pop_chanedit(chans, tmpcom{:}); % apply modification to channel structure
                    if iscell(tmpcom{2}) & (length(tmpcom{2}) == 2) & isstr(tmpcom{2}{2}) & strcmpi(tmpcom{2}{2}, 'gui'),
                        tmpcom = { tmpcom{1} tmpcom{2}{1} };
                    end;
                    if iscell(tmpcom{2}) & (length(tmpcom{2}) > 0) & isstr(tmpcom{2}{1}) & strcmpi(tmpcom{2}{1}, 'chancenter')
                        tmpcom = { tmpcom{1} { tmpcom{2}{1} tmpcom{2}{2} 0 } };
                    end;
                    totaluserdat = { totaluserdat{:} tmpcom{:} };
                catch
                    errordlg2(lasterr, 'Channel location error');
                    returnmode = 'no';
                end;	
                evalin('base', 'clear comtmp');
            end;
            chans = checkchans(chans, allfields);
		end;
		
		% handle arguments
		% ----------------
		if strcmp(returnmode, 'retuninginputui')
			ingui = 0;
			if nbchan ~= 0 & nbchan ~= length(chans)
				if ~popask(strvcat(['The number of channel (' int2str(length(chans)) ') does not correspond to the'], ...
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
            if ~isequal(tmpshrink, shrinkfact)
				if  isstr(tmpshrink)  chans(1).shrink = num2str(tmpshrink);
				else                  chans(1).shrink = tmpshrink;
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
	evalin('base', 'clear global chantmp; clear chantmp;');
else 
     args = varargin;
	 % no interactive inputs
	 % scan all the fields of g
	 % ------------------------
	 for curfield = 1:2:length(args)
		 switch lower(args{curfield})
          case 'forcelocs',
           if ~isempty(args{curfield+1})
               chans = forcelocs(chans,args{curfield+1});
               disp('Convert XYZ coordinates to spherical and polar');
           end;
          case 'eval',
           chantmp = chans;
           eval( args{curfield+1} );
           chans = chantmp;
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
               if ~strcmpi(tmpButtonName, 'Yes'), return; end;
		   end;
		   switch method
			case 'chancenter',
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
             disp('Note: automatically convert XYZ coordinates to spherical and polar');
             chans = convertlocs(chans, 'cart2all');
            otherwise
             chans = convertlocs(chans, method, 'verbose', 'on');
		   end;
		  case 'transform'
		   try, tmpoper = args{curfield+1}; catch, return; end;
		   if isempty(deblank(tmpoper)), return; end;
		   if iscell(tmpoper), tmpoper = tmpoper{1}; end;
           tmpoper = [ tmpoper ';' ];
		   if isempty(findstr(tmpoper, 'chans'))
			   try, X          = cell2mat({chans.X});          catch, X(1:length(chans)) = NaN; end;
			   try, Y          = cell2mat({chans.Y});          catch, Y(1:length(chans)) = NaN; end;
			   try, Z          = cell2mat({chans.Z});          catch, Z(1:length(chans)) = NaN; end;
			   try, theta      = cell2mat({chans.theta});      catch, theta(1:length(chans)) = NaN; end;
			   try, radius     = cell2mat({chans.radius});     catch, radius(1:length(chans)) = NaN; end;
			   try, sph_theta  = cell2mat({chans.sph_theta});  catch, sph_theta(1:length(chans)) = NaN; end;
			   try, sph_phi    = cell2mat({chans.sph_phi});    catch, sph_phi(1:length(chans)) = NaN; end;
			   try, sph_radius = cell2mat({chans.sph_radius}); catch, sph_radius(1:length(chans)) = NaN; end;
			   eval(tmpoper)
			   chans = struct('labels', { chans.labels }, 'X', mat2cell(X), 'Y', mat2cell(Y), 'Z', mat2cell(Z), ...
							  'theta', mat2cell(theta), 'radius', mat2cell(radius), 'sph_theta', mat2cell(sph_theta), ...
							  'sph_phi', mat2cell(sph_phi), 'sph_radius', mat2cell(sph_radius));
               if     ~isempty(findstr(tmpoper, 'X')),          chans = convertlocs(chans, 'cart2all'); end;
               if     ~isempty(findstr(tmpoper, 'Y')),          chans = convertlocs(chans, 'cart2all'); end;
               if     ~isempty(findstr(tmpoper, 'Z')),          chans = convertlocs(chans, 'cart2all'); end;
               if     ~isempty(findstr(tmpoper, 'sph_theta')),  chans = convertlocs(chans, 'sph2all');
               elseif ~isempty(findstr(tmpoper, 'theta')),      chans = convertlocs(chans, 'topo2all'); end;
               if     ~isempty(findstr(tmpoper, 'sph_phi')),    chans = convertlocs(chans, 'sph2all');  end;
               if     ~isempty(findstr(tmpoper, 'sph_radius')), chans = convertlocs(chans, 'sph2all');
               elseif ~isempty(findstr(tmpoper, 'radius')),     chans = convertlocs(chans, 'topo2all'); end;
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
           
		  case { 'insert' 'add' }
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
           
		  case 'append'
		   tmpargs = args{ curfield+1 };
		   allfields = fieldnames(chans);
		   if length( tmpargs ) < length(allfields)+1
			   error('pop_chanedit: not enough arguments to change all field values');
		   end;
		   num = tmpargs{1};
		   chans(end+1) = chans(end);
		   chans(num+2:end) = chans(num+1:end-1);
		   for index = 1:length( allfields )
               eval([ 'chans(' int2str(num+1) ').' allfields{index} '=' reformat(tmpargs{index+1}) ';' ]);
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
                   chans = readlocs(tmpargs{:});
			   end;
		   end;
           
		  case 'lookup'
           tmplocs = readlocs('Standard-10-20-Cap81.ced', 'filetype', 'chanedit');
           [tmp ind1 ind2] = intersect(lower({ tmplocs.labels }), lower({ chans.labels }));
           if ~isempty(tmp)
               [ind2 ind3] = sort(ind2);
               ind1 = ind1(ind3);
               for index = 1:length(ind2)
                   chans(ind2(index)).theta      = tmplocs(ind1(index)).theta;
                   chans(ind2(index)).radius     = tmplocs(ind1(index)).radius;
                   chans(ind2(index)).X          = tmplocs(ind1(index)).X;
                   chans(ind2(index)).Y          = tmplocs(ind1(index)).Y;
                   chans(ind2(index)).Z          = tmplocs(ind1(index)).Z;
                   chans(ind2(index)).sph_theta  = tmplocs(ind1(index)).sph_theta;
                   chans(ind2(index)).sph_phi    = tmplocs(ind1(index)).sph_phi;
                   chans(ind2(index)).sph_radius = tmplocs(ind1(index)).sph_radius;
               end;
               tmpdiff = setdiff([1:length(chans)], ind2);
               if ~isempty(tmpdiff)
                   fprintf('Channel lookup: no location for ');
                   for index = 1:(length(tmpdiff)-1)
                       fprintf('%s,', chans(tmpdiff(index)).labels);
                   end;
                   fprintf('%s\nSend us standard location for your channels at eeglab@sccn.ucsd.edu\n', ...
                           chans(tmpdiff(end)).labels);
               end;
           end;
          case 'eval'
		   tmpargs = args{ curfield+1 };
           eval(tmpargs);
		  
          case 'save'
		   tmpargs = args{ curfield+1 };
		   if isempty(tmpargs), return; end;
		   fid = fopen(tmpargs, 'w');
		   if fid ==-1, error('Cannot open file'); end;
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

if isfield(chans, 'sph_phi_besa'  ), chans = rmfield(chans, 'sph_phi_besa'); end;
if isfield(chans, 'sph_theta_besa'), chans = rmfield(chans, 'sph_theta_besa'); end;
return;

function str = commandwarning(str);

% format the output field
% -----------------------
function strval = reformat( val )
    if isnumeric(val) & isempty(val), val = '[]'; end;
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
function [chans, shrinkfact]= checkchans(chans, fields);
	shrinkfact = 'off';
	if isfield(chans, 'shrink'), 
		shrinkfact = chans(1).shrink; 
		chans = rmfield(chans, 'shrink');
        if isstr(shrinkfact) & ~isempty(str2num(shrinkfact)), shrinkfact = str2num(shrinkfact); end;
	end;
    for index = 1:length(fields)
        if ~isfield(chans, fields{index})
            chans = setfield(chans, {1}, fields{index}, []);
        end;
    end;

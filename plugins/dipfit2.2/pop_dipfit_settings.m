% pop_dipfit_settings() - select global settings for dipole fitting through a pop up window
%
% Usage:
%   >> OUTEEG = pop_dipfit_settings ( INEEG ); % pop up window
%   >> OUTEEG = pop_dipfit_settings ( INEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Inputs:
%   INEEG	input dataset
%
% Optional inputs:
%   'hdmfile'  - [string] file containing a head model compatible with
%                the Fieldtrip dipolefitting() function ("vol" entry)
%   'mrifile'  - [string] file containing an anatomical MR head image. 
%                The MRI must be normalized to the MNI brain. See the .mat 
%                files used by the sphere and boundary element models
%                (For instance, select the sphere model and study 'EEG.dipfit'). 
%                If SPM2 software is installed, dipfit will be able to read 
%                most MRI file formats for plotting purposes (.mnc files, etc...). 
%                To plot dipoles in a subject MRI, first normalize the MRI 
%                to the MNI brain using SPM2.
%   'coordformat' - ['MNI'|'Spherical'] Coordinates returned by the selected
%                head model. May be MNI coordinates or spherical coordinates
%                (For spherical coordinates, the head radius is assumed to be 85 mm.
%   'chanfile' - [string] template channel locations file. (This function will
%                check whether your channel locations file is compatible with 
%                your selected head model).
%   'chansel'  - [int. vector] indices of channels to use for dipole fitting. 
%                {default: all} ????????????????
%   'coord_transform' - [float array] Talairach transformation matrix for
%                       aligning the dataset channel locations to the selected 
%                       head model.
%   'electrodes'      - [integer array] indices of channels to include
%                       in the dipole model. {default: all}
% Outputs:
%   OUTEEG	output dataset
%
% Author: Arnaud Delorme, SCCN, La Jolla 2003-
%         Robert Oostenveld, SMI/FCDC, Nijmegen 2003

% MEG flag:
%   'gradfile' - [string] file containing gradiometer locations
%                ("gradfile" parameter in Fieldtrip dipolefitting() function)

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl

% Copyright (C) 2003 arno@salk.edu, Arnaud Delorme, SCCN, La Jolla 2003-2005
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
% Revision 1.26  2006/09/19 19:52:33  toby
% further text and help edits
%
% Revision 1.25  2006/09/19 19:41:05  toby
% further text edits
%
% Revision 1.24  2006/09/19 19:32:45  toby
% edited gui text -sm & tf
%
% Revision 1.23  2006/05/17 18:32:51  arno
% fixing history
%
% Revision 1.22  2006/04/13 17:33:02  arno
% coreg checkbox disable
%
% Revision 1.21  2006/04/13 17:30:48  arno
% same
%
% Revision 1.20  2006/04/13 17:30:09  arno
% same
%
% Revision 1.19  2006/04/13 17:28:43  arno
% default values for model
%
% Revision 1.18  2006/03/07 00:46:05  arno
% final edits
%
% Revision 1.17  2006/03/07 00:44:29  arno
% same
%
% Revision 1.16  2006/03/07 00:42:14  arno
% same
%
% Revision 1.15  2006/03/07 00:41:10  arno
% update help message
%
% Revision 1.14  2006/03/07 00:37:47  arno
% disable edit boxes
%
% Revision 1.13  2006/01/19 22:06:43  arno
% correcting header
%
% Revision 1.12  2006/01/12 23:07:32  arno
% now processing chaninfo
%
% Revision 1.11  2006/01/11 00:13:38  arno
% adding help message for coregister
%
% Revision 1.10  2006/01/10 00:44:33  arno
% allowing coregistration
%
% Revision 1.9  2005/03/21 19:03:00  arno
% adding one line
%
% Revision 1.8  2005/03/17 17:33:41  arno
% warning if dipole info already present
%
% Revision 1.7  2005/03/16 02:32:39  arno
% detect channels with no coordinates
%
% Revision 1.6  2005/03/11 16:07:07  arno
% header
%
% Revision 1.5  2005/03/10 19:34:58  arno
% backward compatibility with old option electrodes
%
% Revision 1.4  2005/03/10 18:48:12  arno
% adding new entry to select coordinate format
%
% Revision 1.3  2005/03/10 18:07:22  arno
% lowercase
%
% Revision 1.2  2005/03/10 18:05:27  arno
% renaming BESA file
%
% Revision 1.1  2005/03/10 18:04:39  arno
% Initial revision
%
% Revision 1.27  2004/01/07 17:02:48  scott
% See List -> List
%
% Revision 1.26  2004/01/07 17:01:32  scott
% ... -> See List
%
% Revision 1.25  2004/01/07 17:00:21  scott
% added to Note
%
% Revision 1.24  2003/12/04 18:25:43  arno
% shell conductance unit
%
% Revision 1.23  2003/10/31 18:06:07  arno
% more edits
%
% Revision 1.22  2003/10/31 17:54:00  arno
% select channels -> omit channels
%
% Revision 1.21  2003/10/30 02:14:51  arno
% gui typo
%
% Revision 1.20  2003/10/29 22:43:37  arno
% wording
%
% Revision 1.19  2003/10/29 02:36:38  arno
% removing 2 lines of GUI
%
% Revision 1.18  2003/10/15 14:47:34  roberto
% removed urchanlocs from electrode projection part
%
% Revision 1.17  2003/10/14 15:56:38  roberto
% before projecting electrodes towards the skin, store the original in urchanlocs
%
% Revision 1.16  2003/08/08 16:57:10  arno
% removing normsphere
%
% Revision 1.15  2003/08/04 22:10:13  arno
% adding warning backtrace
%
% Revision 1.14  2003/08/04 22:03:32  arno
% adding normsphere option
%
% Revision 1.13  2003/08/01 13:50:51  roberto
% changed the gui, implemented fit-electrodes-to-sphere, modifications in vol.r/c/o are accepter
%
% Revision 1.12  2003/07/01 22:11:59  arno
% removing debug message
%
% Revision 1.11  2003/06/30 02:11:37  arno
% debug argument check
%
% Revision 1.10  2003/06/30 01:20:58  arno
% 4sphere -> 4spheres
%
% Revision 1.9  2003/06/30 01:18:46  arno
% copying new version
%
% Revision 1.5  2003/06/16 15:32:51  arno
% reprograming interface, programing history
%
% Revision 1.4  2003/03/12 10:32:50  roberto
% added 4-sphere volume model similar to BESA
%
% Revision 1.3  2003/03/06 15:58:28  roberto
% *** empty log message ***
%
% Revision 1.1  2003/02/24 10:06:08  roberto
% Initial revision
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUTEEG, com] = pop_dipfit_settings ( EEG, varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help pop_dipfit_settings;
   return;
end;

OUTEEG = EEG;
com = '';

% get the default values and filenames
dipfitdefs;

if nargin < 2
    
    % detect DIPFIT1.0x structure
    % ---------------------------
    if isfield(EEG.dipfit, 'vol')
        str = [ 'Dipole information structure from DIPFIT v1.02 detected.' ...
                'Keep or erase the old dipole information including dipole locations? ' ...
                'In either case, a new dipole model can be constructed.' ];
        
        tmpButtonName=questdlg2( strmultiline(str, 60), 'Old DIPFIT structure', 'Keep', 'Erase', 'Keep');
        if strcmpi(tmpButtonName, 'Keep'), return; end;       

    elseif isfield(EEG.dipfit, 'hdmfile')
        % detect previous DIPFIT structure
        % --------------------------------
        str = [ 'Dipole information and settings are present in the dataset. ' ...
                'Keep or erase this information?' ];
        tmpButtonName=questdlg2( strmultiline(str, 60), 'Old DIPFIT structure', 'Keep', 'Erase', 'Keep');
        if strcmpi(tmpButtonName, 'Keep'), return; end;       
    end;    
    
    % define the callbacks for the buttons
    % -------------------------------------
    cb_selectelectrodes = [ 'tmp = select_channel_list({EEG.chanlocs.label}, ' ...
                            'eval(get(findobj(gcbf, ''tag'', ''elec''), ''string'')));' ...
                            'set(findobj(gcbf, ''tag'', ''elec''), ''string'',[''[''  num2str(tmp) '']''])' ]; % did not work
    cb_selectelectrodes = 'set(findobj(gcbf, ''tag'', ''elec''), ''string'', int2str(pop_chansel({EEG.chanlocs.labels})));';
    cb_volmodel = [ 'tmpdat = get(gcbf, ''userdata'');' ... 
                    'tmpind = get(gcbo, ''value'');' ... 
                    'set(findobj(gcbf, ''tag'', ''radii''),   ''string'', num2str(tmpdat{tmpind}.r,3));' ...
                    'set(findobj(gcbf, ''tag'', ''conduct''), ''string'', num2str(tmpdat{tmpind}.c,3));' ...
                    'clear tmpdat tmpind;' ];
    cb_changeradii   = [  'tmpdat = get(gcbf, ''userdata'');' ...
                          'tmpdat.vol.r = str2num(get(gcbo, ''string''));' ...
                          'set(gcf, ''userdata'', tmpdat)' ];
    cb_changeconduct = [  'tmpdat = get(gcbf, ''userdata'');' ...
                          'tmpdat.vol.c = str2num(get(gcbo, ''string''));' ...
                          'set(gcf, ''userdata'', tmpdat)' ];
    cb_changeorigin  = [  'tmpdat = get(gcbf, ''userdata'');' ...
                          'tmpdat.vol.o = str2num(get(gcbo, ''string''));' ...
                          'set(gcf, ''userdata'', tmpdat)' ];
    % cb_fitelec = [ 'if get(gcbo, ''value''),' ...
    %                '  set(findobj(gcbf, ''tag'', ''origin''), ''enable'', ''off'');' ...
    %                'else' ...
    %                '  set(findobj(gcbf, ''tag'', ''origin''), ''enable'', ''on'');' ...
    %                'end;' ];
    valmodel     = 1;
    nocoregvalue = 0;
    coregvals    = '';
    if isfield(EEG.chaninfo, 'filename')
        if ~isempty(findstr(lower(EEG.chaninfo.filename), 'standard-10-5-cap385')), nocoregvalue = 1; end;
        if ~isempty(findstr(lower(EEG.chaninfo.filename), 'standard_1005')),        coregvals = num2str(template_models{valmodel}{5}); valmodel = 2; end;
    end;
        
    userdata    = [];
    
    geomvert = [2 1 1 1 1 1 1 1 1 1 1];
    
    geomhorz = {
        [1 2] 
        [1]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.9 0.1 ]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.5 0.5 ]
        [1 1.3 0.5 0.5 ]
        [1]
        [1]
        [1] };
    
    % define each individual graphical user element
    comhelp1 = [ 'warndlg2(strvcat(''The two default head models are in ''standard_BEM'' and ''standard_BESA'''',' ...
                 ''' sub-folders in the DIPFIT2 plugin folder, and may be modified there.''), ''Model type'');' ];
    comhelp3 = [ 'warndlg2(strvcat(''Any MR image normalized to the MNI brain model may be used for plotting'',' ...
                 '''(see the DIPFIT 2.0 tutorial for more information)''), ''Model type'');' ];
    comhelp2 = [ 'warndlg2(strvcat(''The template location file associated with the head model'',' ...
                 '''you are using must be entered (see tutorial).''), ''Template location file'');' ];
    commandload1 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''model''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    commandload2 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''meg''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    commandload3 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''mri''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];    
    cb_selectcoreg = [ 'tmpmodel = get( findobj(gcbf, ''tag'', ''model''), ''string'');' ...
                       'tmploc2  = get( findobj(gcbf, ''tag'', ''meg'')  , ''string'');' ...
                       'tmploc1  = get( gcbo, ''userdata'');' ...
                       'tmptransf = get( findobj(gcbf, ''tag'', ''coregtext''), ''string'');' ...
                       '[tmp tmptransf] = coregister(tmploc1{1}, tmploc2, ''mesh'', tmpmodel,' ...
                       '                       ''transform'', str2num(tmptransf), ''chaninfo1'', tmploc1{2}, ''helpmsg'', ''on'');' ...
                       'if ~isempty(tmptransf), set( findobj(gcbf, ''tag'', ''coregtext''), ''string'', num2str(tmptransf)); end;' ...
                       'clear tmpmodel tmploc2 tmploc1 tmp tmptransf;' ];
    
    dipfitdefs; % contains template_model
    setmodel = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                 'tmpval = get(gcbo, ''value'');' ...
                 'set(findobj(gcbf, ''tag'', ''model''), ''string'', tmpdat{tmpval}{1});' ...
                 'set(findobj(gcbf, ''tag'', ''coord''), ''value'' , fastif(strcmpi(tmpdat{tmpval}{2},''MNI''),2,1));' ...
                 'set(findobj(gcbf, ''tag'', ''mri''  ), ''string'', tmpdat{tmpval}{3});' ...
                 'set(findobj(gcbf, ''tag'', ''meg''), ''string'', tmpdat{tmpval}{4});' ...
                 'set(findobj(gcbf, ''tag'', ''coregcheckbox''), ''value'', 0);' ...
                 'if tmpval < 3,' ...
                 '  set(findobj(gcbf, ''userdata'', ''editable''), ''enable'', ''off'');' ...
                 'else,' ...
                 '  set(findobj(gcbf, ''userdata'', ''editable''), ''enable'', ''on'');' ...
                 'end;' ];
        
    elements  = { ...
        { 'style' 'text'        'string'  [ 'Head model (click to select)' 10 '' ] } ...
        { 'style' 'listbox'     'string'  'Spherical Four-Shell (BESA) | Boundary Element Model (MNI) | Custom Head' ... 
                                'callback' setmodel 'value' valmodel } { } ...
        { 'style' 'text'        'string' 'Head model file' } ...
        { 'style' 'edit'        'string' template_models{valmodel}{1}  'tag'      'model' 'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'pushbutton'  'string' 'Browse'    'callback' commandload1       'userdata' 'editable' 'enable' 'off' } ...
        { 'style' 'pushbutton'  'string' 'Help'      'callback' comhelp1 } ...
        { 'style' 'text'        'string' 'Ourput coordinates' } ...
        { 'style' 'listbox'     'string' 'spherical (head radius 85 mm)|MNI' 'tag' 'coord' ...
          'value' fastif(strcmpi(template_models{valmodel}{2},'MNI'),2,1)  'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'text'        'string' 'Click to select' } { } ...
        { 'style' 'text'        'string' 'MRI' } ...
        { 'style' 'edit'        'string' template_models{valmodel}{3} 'tag'      'mri' } ...
        { 'style' 'pushbutton'  'string' 'Browse'       'callback' commandload3 } ...
        { 'style' 'pushbutton'  'string' 'Help'         'callback' comhelp3 } ...
        { 'style' 'text'        'string' 'Template channel locations file' } ...
        { 'style' 'edit'        'string' template_models{valmodel}{4} 'tag'      'meg'  'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'pushbutton'  'string' 'Browse'       'callback' commandload2  'userdata' 'editable' 'enable' 'off'} ...
        { 'style' 'pushbutton'  'string' 'Help'         'callback' comhelp2 } ...
        { 'style' 'text'        'string' 'Co-register chan. locs. with head model' } ...
        { 'style' 'edit'        'string' coregvals 'tag' 'coregtext' } ...
        { 'style' 'pushbutton'  'string' 'Manual Co-Reg.' 'callback' cb_selectcoreg 'userdata' { EEG.chanlocs EEG.chaninfo } } ... 
        { 'style' 'checkbox'    'string' 'No Co-Reg.'    'tag' 'coregcheckbox' 'value' nocoregvalue } ... 
        { 'style' 'text'        'string' 'Channels to omit from dipole fitting' } ...
        { 'style' 'edit'        'string' ''             'tag' 'elec' } ...
        { 'style' 'pushbutton'  'string' 'List' 'callback' cb_selectelectrodes } { } ...
        { } ...
        { 'style' 'text'        'string' 'Note: Check that the channel locations are on the surface of the head model' } ...
        { 'style' 'text'        'string' '(To do this: ''Set head radius'' to about 85 in the channel editor).' } ...
                };
    
    % plot GUI and protect parameters
    % -------------------------------
    optiongui = { 'geometry', geomhorz, 'uilist', elements, 'helpcom', 'pophelp(''pop_dipfit_settings'')', ...
                  'title', 'Dipole fit settings - pop_dipfit_settings()', ...
                  'userdata', template_models, 'geomvert', geomvert };
	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', 'noclose', optiongui{:});
    if isempty(result), return; end;
    if ~isempty(get(0, 'currentfigure')) currentfig = gcf; else return; end;
    
    while test_wrong_parameters(currentfig)
    	[result, userdat2, strhalt, outstruct] = inputgui( 'mode', currentfig, optiongui{:});
        if isempty(result), return; end;
    end;
    close(currentfig);

    % decode GUI inputs
    % -----------------
    options = {};
    options = { options{:} 'hdmfile'      result{2} };
    options = { options{:} 'coordformat'  fastif(result{3} == 2, 'MNI', 'Spherical') };
    options = { options{:} 'mrifile'      result{4} };
    options = { options{:} 'chanfile'     result{5} };
    if ~result{7}, options = { options{:} 'coord_transform' str2num(result{6}) }; end;
    options = { options{:} 'chansel'      setdiff(1:EEG.nbchan, str2num(result{8})) };

else
    options = varargin;
end

g = finputcheck(options, { 'hdmfile'  'string'    []         '';
                                 'mrifile'  'string'    []         '';
                                 'chanfile' 'string'    []         '';
                                 'chansel'  'integer'   []         [1:EEG.nbchan];
                                 'electrodes' 'integer'   []         [];
                                 'coord_transform' 'real' []         [];
                                 'coordformat' 'string'    { 'MNI' 'spherical' } 'MNI' });
if isstr(g), error(g); end;

OUTEEG = rmfield(OUTEEG, 'dipfit');
OUTEEG.dipfit.hdmfile     = g.hdmfile;
OUTEEG.dipfit.mrifile     = g.mrifile;
OUTEEG.dipfit.chanfile    = g.chanfile;
OUTEEG.dipfit.chansel     = g.chansel;
OUTEEG.dipfit.coordformat     = g.coordformat;
OUTEEG.dipfit.coord_transform = g.coord_transform;
if ~isempty(g.electrodes), OUTEEG.dipfit.chansel = g.electrodes; end;

% removing channels with no coordinates
% -------------------------------------
[tmpeloc labels Th Rd indices] = readlocs(EEG.chanlocs);
if length(indices) < length(EEG.chanlocs)
    disp('Warning: Channels removed from dipole fitting no longer have location coordinates!');
    OUTEEG.dipfit.chansel = intersect( OUTEEG.dipfit.chansel, indices);
end;

% checking electrode configuration
% --------------------------------
if 0
    disp('Checking the electrode configuration');
    tmpchan          = readlocs(OUTEEG.dipfit.chanfile);
    [tmp1 ind1 ind2] = intersect( lower({ tmpchan.labels }), lower({ OUTEEG.chanlocs.labels }));
    if isempty(tmp1)
        disp('No channel labels in common found between template and dataset channels');
        if ~isempty(findstr(OUTEEG.dipfit.hdmfile, 'BESA'))
            disp('Use the channel editor to fit a head sphere to your channel locations.');
            disp('Check for inconsistency in dipole info.');
        else
            disp('Results using standard BEM model are INACCURATE when the chan locations are not on the head surface!');
        end;
    else % common channels: performing best transformation
        TMP = OUTEEG;
        elec1 = eeglab2fieldtrip(TMP, 'elec');
        elec1 = elec1.elec;
        TMP.chanlocs = tmpchan;
        elec2 = eeglab2fieldtrip(TMP, 'elec');
        elec2 = elec2.elec;
        cfg.elec     = elec1;
        cfg.template = elec2;
        cfg.method   = 'warp';
        elec3 = electrodenormalize(cfg);
    
        % convert back to EEGLAB format
        OUTEEG.chanlocs = struct( 'labels', elec3.label, ...
                              'X'     , mat2cell(elec3.pnt(:,1)'), ...
                              'Y'     , mat2cell(elec3.pnt(:,2)'), ...
                              'Z'     , mat2cell(elec3.pnt(:,3)') );
        OUTEEG.chanlocs = convertlocs(OUTEEG.chanlocs, 'cart2all');
    end;
    
end;

com = sprintf('%s = pop_dipfit_settings( %s, %s);', inputname(1), inputname(1), vararg2str(options));

% test for wrong parameters
% -------------------------
function bool = test_wrong_parameters(hdl)

    coreg1 = get( findobj( hdl, 'tag', 'coregtext')    , 'string' );
    coreg2 = get( findobj( hdl, 'tag', 'coregcheckbox'), 'value' );
    
    bool = 0;
    if coreg2 == 0 & isempty(coreg1)
         bool = 1; warndlg2(strvcat('You must co-register your channel locations', ...
                                    'with the head model (Press buttun, "Manual Co-Reg".', ...
                                    'and follow instructions); To bypass co-registration,', ...
                                    'check the checkbox " No Co-Reg".'), 'Error');
    end;

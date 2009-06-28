% pop_chanedit() - Edit the channel locations structure of an EEGLAB dataset,
%                  EEG.chanlocs. For structure location and file formats,
%                  see >> help readlocs
%
% Usage:    >> newchans = pop_chanedit( EEG, 'key1', value1, ...
%                        'key2', value2, ... ); % edit dataset containing chanlocs
%           >> [ newchans options ] = pop_chanedit( chanlocs, 'key1', value1, ...
%                        'key2', value2, ... ); % edit separate chanlocs
%                        struct
% Graphic interface:
%   "Channel information ('field name')" - [edit boxes] display channel field
%                   contents for the current channel. Command line equivalent
%                   to modify these fields: 'transform'
%   "Opt. 3D center" - [button] optimally re-center 3-D channel coordinates. Uses
%                 chancenter(). Command line equivalent: 'convert', { 'chancenter'
%                 [xc yc zc] }, [xc yc zc] being the center of the sphere. Use []
%                 to find the center of the best fitting sphere.
%   "Rotate axis" - [button] force one electrode to one position and rotate the other
%                 electrodes accordingly. Command line equivalent: 'forcelocs'.
%   "Transform axis" - [button] perform any operation on channel fields. Command
%                 line equivalent: 'transform'.
%   "Xyz->polar & sph." - [button] convert 3-D cartesian coordinates to polar and
%                 3-D spherical coordinates. This is useful when you edit the
%                 coordinates manually. Command line equivalent: 'convert', 'cart2all'.
%   "Sph.->polar & xyz" - [button] convert 3-D spherical coordinates to polar and
%                 3-D cartesian coordinates. Command line equivalent: 'convert', 'sph2all'.
%   "Polar->sph & xyz" - [button] convert 2-D polar coordinates to 3-D spherical and
%                 3-D cartesian coordinates. Command line equivalent: 'convert', 'topo2all'.
%                 Note that if spherical radii are absent, they are forced to 1.
%   "Set head radius" - [button] change head size radius. This is useful
%                 to make channels location compatible with a specified spherical model.
%                 Command line equivalent: 'headrad'.
%   "Set channel types" - [button] set channel type names for a range of data channels.
%   "Delete chan" - [button] delete channel. Command line equivalent: 'delete'.
%   "Insert chan" - [button] insert channel before current channel.
%                 Command line equivalent: 'insert'.
%   "<<" - [button] scroll channel backward by 10.
%   "<" - [button] scroll channel backward by 1.
%   ">" - [button] scroll channel forward by 1.
%   ">>" - [button] scroll channel forward by 10.
%   "Append chan" - [button] append channel after the current channel.
%                 Command line equivalent: 'append'.
%   "Plot 2D"     - [button] plot channel locations in 2-D using topoplot()
%   "Plot radius [value (0.2-1.0), []=auto)" - [edit box] default plotting radius
%                 in 2-D polar views. This does NOT affect channel locations; it
%                 is only used for visualization. This parameter is attached to the
%                 chanlocs structure and is then used in all 2-D scalp topoplots.
%                 Default -> to data limits. Command line equivalent: 'plotrad'.
%   "Nose along +X" - [list] Indicate the direction of the nose. This information
%                 is used in functions like topoplot(), headplot() and dipplot().
%                 Command line equivalent: 'nosedir'.
%   "Plot 3D"     - [button] plot channel positions in 3-D using plotchans3d()
%   "Read locations" - [button] read location file using readlocs()
%                 Command line equivalent: 'load'.
%   "Read help"   - [button] display readlocs() function help.
%   "Save .ced"   - [button] save channel locations in native EEGLAB ".ced" format.
%                 Command line equivalent: 'save'.
%   "Save others" - [button] save channel locations in other formats using
%                 pop_writelocs() (see readlocs() for available channel formats).
%   "Cancel"      - [button] cancel all editing.
%   "Help"        - [button] display this help message.
%   "OK"          - [button] save edits and propagate to parent.
%
% Inputs:
%   EEG      - EEG dataset
%   chanlocs - EEG.chanlocs structure
%
% Optional inputs:
%   'convert'     - {conversion_type [args]} Conversion type may be: 'cart2topo'
%                   'sph2topo', 'topo2sph', 'sph2cart', 'cart2sph', or 'chancenter'.
%                   See help messages for these functions. Args are only relevant
%                   for 'chancenter'. More info is given in the graphic interface
%                   description above.
%   'transform'   - String command for manipulating arrays. 'chan' is full channel
%                   info. Fields that can be manipulated are 'labels', 'theta'
%                   'radius' (polar angle and radius), 'X', 'Y', 'Z' (cartesian
%                   3-D) or 'sph_theta', 'sph_phi', 'sph_radius' for spherical
%                   horizontal angle, azimuth and radius.
%                   Ex: 'chans(3) = chans(14)', 'X = -X' or a multi-step transform
%                   with steps separated by ';': Ex. 'TMP = X; X = Y; Y = TMP'
%   'changechan'  - {number value1 value2 value3 ...} Change the values of all fields
%                   for the given channel number, mimimally {num label theta radius}.
%                   Ex: 'changechan' {12 'PXz' -90 0.30}
%   'changefield' - {number field value} Change field value for channel number number.
%                   Ex: {34 'theta' 320.4}.
%   'insert'      - {number label theta radius X Y Z sph_theta sph_phi sph_radius }
%                   Insert new channel and specified values before the current channel
%                    number. If the number of values is less than 10, remaining
%                   fields will be 0. (Previously, this parameter was termed 'add').
%   'append'      - {num label theta radius X Y Z sph_theta sph_phi sph_radius }
%                   same as 'insert' (above) but insert the the new channel after
%                   the current channel number.
%   'delete'      - [int_vector] Vector of channel numbers to delete.
%   'forcelocs'   - [cell] call forcelocs() to force a particular channel to be at a
%                   particular location on the head sphere; rotate other channels
%                   accordingly.
%   'skirt'       - Topographical polar skirt factor (see >> help topoplot)
%   'shrink'      - Topographical polar shrink factor (see >> help topoplot)
%   'load'        - [filename|{filename, 'key', 'val'}] Load channel location file
%                   optional arguments (such as file format) to the function
%                   readlocs() can be specified if the input is a cell array.
%   'save'        - 'filename' Save text file with channel info.
%   'eval'        - [string] evaluate string ('chantmp' is the name of the channel
%                   location structure).
%   'headrad'     - [float] change head radius.
%   'lookup'      - [string] look-up channel numbers for standard locations in the
%                   channel location file given as input.
%
% Outputs:
%   newchans      - new EEGLAB channel locations structure or EEG dataset with updated
%                   channel location structures EEG.chanlocs, EEG.urchanlocs, EEG.chaninfo
%   options       - structure containing plotting options (equivalent to EEG.chaninfo)
%
% Ex:    EEG = pop_chanedit(EEG,'load', { 'dummy.elp' 'elp' }, 'delete', [3 4], ...
%                       'convert', { 'xyz->polar' [] -1 1 }, 'save', 'mychans.loc' )
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
% Revision 1.177  2009/06/28 05:49:56  arno
% Adding reference and reprogramming pop_chanedit
%
% Revision 1.176  2008/11/11 02:52:26  arno
% Do not test if structure any more
%
% Revision 1.175  2008/02/15 16:32:20  arno
% adding more recognized channel types
%
% Revision 1.174  2007/09/11 15:03:05  arno
% fix channel removal (fiducials etc...)
%
% Revision 1.173  2007/08/24 01:45:48  arno
% remove shiftdatachan obsolete parameter
%
% Revision 1.172  2007/08/24 01:31:00  scott
% many help msg edits -- Arno, please find and fix ????? if necessary
%
% Revision 1.171  2007/08/24 00:53:36  arno
% dealing with extra channels
%
% Revision 1.170  2007/08/23 23:37:23  arno
% handle reference
%
% Revision 1.169  2007/08/23 19:25:40  arno
% debug if no argument
%
% Revision 1.168  2007/08/16 19:02:23  arno
% different template_models format
%
% Revision 1.167  2007/08/15 19:26:13  arno
% same
%
% Revision 1.166  2007/08/15 19:25:39  arno
% file filter
%
% Revision 1.165  2006/12/06 23:39:01  arno
% fiducial fix
%
% Revision 1.164  2006/12/05 23:26:42  arno
% take into account fiducials from command line
%
% Revision 1.163  2006/12/05 21:40:53  arno
% fixing warning message nbchan difference...
%
% Revision 1.162  2006/11/22 18:50:02  arno
% fix chaninfo update from command line
%
% Revision 1.161  2006/11/21 22:43:20  arno
% debug lookup for BEM file
%
% Revision 1.160  2006/11/10 01:58:52  arno
% text
%
% Revision 1.159  2006/11/06 21:45:34  arno
% fix .nodatchans for fieldtrip compatibility
%
% Revision 1.158  2006/11/03 23:32:01  arno
% fix lookup gui, text etc...
%
% Revision 1.157  2006/11/03 23:22:17  arno
% same
%
% Revision 1.156  2006/11/03 23:20:50  arno
% changing message fro wrong number of channels
%
% Revision 1.155  2006/11/03 22:05:16  arno
% update message
%
% Revision 1.154  2006/10/25 21:25:07  arno
% fixing integer channel type
%
% Revision 1.153  2006/10/04 16:29:52  arno
% fix nosedir  history
%
% Revision 1.152  2006/10/04 14:40:44  arno
% fixed integer channel type
%
% Revision 1.150  2006/04/17 15:47:54  scott
% help message formatting
%
% Revision 1.149  2006/04/10 08:30:33  scott
% dos2unix
%
% Revision 1.148  2006/04/09 02:37:06  scott
% window text detail
%
% Revision 1.147  2006/01/26 21:04:48  arno
% nothing
%
% Revision 1.146  2006/01/21 00:13:41  arno
% saving filename also when lookup
%
% Revision 1.145  2006/01/19 19:44:34  arno
% fix popupmenu
%
% Revision 1.144  2006/01/19 19:40:59  arno
% replace listbox by pop-up menu
%
% Revision 1.143  2006/01/19 18:43:30  arno
% remove datachan
%
% Revision 1.142  2006/01/19 18:40:14  arno
% do not erase channel type when lookup
%
% Revision 1.141  2006/01/19 00:37:58  arno
% fixing default elp BESA
%
% Revision 1.140  2006/01/12 23:44:43  arno
% removing nodatchans field from params after transfering nodatchans to structure
%
% Revision 1.139  2006/01/12 23:01:37  arno
% fixing fiducials
%
% Revision 1.138  2006/01/12 22:55:21  arno
% nothing
% ./
%
% Revision 1.137  2006/01/12 22:37:40  arno
% allowing fiducials to be processed
%
% Revision 1.136  2005/12/03 00:15:34  arno
% fix urchanlocs
%
% Revision 1.135  2005/11/09 23:19:37  arno
% backing file content in info structure
%
% Revision 1.134  2005/11/09 22:44:35  arno
% plotting urchans correctly
%
% Revision 1.133  2005/11/09 00:42:30  arno
% urchanlocs
%
% Revision 1.132  2005/10/10 17:35:46  arno
% fix problem when looking-up channel location file
%
% Revision 1.131  2005/09/27 21:56:25  arno
% change aspect ratio for windows; remove channel locs when looking up location
%
% Revision 1.130  2005/09/13 16:12:09  scott
% help msg
%
% Revision 1.129  2005/09/08 21:40:39  arno
% fix return call for multiple datasets
%
% Revision 1.128  2005/08/04 19:02:29  arno
% fix problem when not returning EEG structure
%
% Revision 1.127  2005/08/02 01:34:32  arno
% fixing chansout
%
% Revision 1.126  2005/07/29 23:37:40  arno
% allowing changing several datasets
%
% Revision 1.125  2005/07/29 23:29:31  arno
% allow processing several structures
%
% Revision 1.124  2005/05/24 17:37:20  arno
% remove cell2mat
%
% Revision 1.123  2005/04/08 23:10:22  arno
% fixing nose orientation when looking BESA file
%
% Revision 1.122  2005/03/30 22:31:17  arno
% debug for Matlab 5.3
%
% Revision 1.121  2005/03/16 02:47:26  arno
% fixing transform when some channels do not have coordinates
%
% Revision 1.120  2005/03/11 18:48:56  arno
% lookup button
%
% Revision 1.119  2005/03/11 01:59:27  arno
% autodetecting BEM file for nose orientation
%
% Revision 1.118  2005/03/09 19:29:42  arno
% updating GUI for channel location file
%
% Revision 1.117  2005/03/08 21:49:26  arno
% deug plotrad
%
% Revision 1.116  2005/03/05 02:20:18  arno
% default shrinkorskirt
%
% Revision 1.115  2005/03/05 02:14:33  arno
% debug chaninfo
%
% Revision 1.114  2005/03/04 23:20:41  arno
% new structure chaninfo etc...
%
% Revision 1.111  2004/11/19 02:07:38  arno
% display warning when deleting channels
%
% Revision 1.110  2004/11/12 21:37:30  arno
% nothing
%
% Revision 1.109  2004/11/10 16:05:15  arno
% fix channel type problem
%
% Revision 1.108  2004/11/05 18:06:31  arno
% putting back orderfields with test
%
% Revision 1.107  2004/10/07 17:28:26  hilit
% cosmetic GUI changes
%
% Revision 1.106  2004/10/07 16:28:02  hilit
% disabled the call to orderfields
%
% Revision 1.105  2004/09/15 01:51:14  arno
% nothing
%
% Revision 1.104  2004/09/14 18:04:00  arno
% implementing GUI to set channel types
%
% Revision 1.103  2004/09/14 17:35:26  arno
% nothing
%
% Revision 1.102  2004/07/07 22:37:14  arno
% debug shrink conversion to plotrad
%
% Revision 1.101  2004/07/07 22:31:22  arno
% shrink message
%
% Revision 1.100  2004/07/07 22:21:13  arno
% debug shrinkorskirt
%
% Revision 1.99  2004/06/09 17:44:55  arno
% lookup coordinate
%
% Revision 1.98  2004/05/22 00:00:00  arno
% nothing
%
% Revision 1.97  2004/05/21 17:54:46  arno
% allowing to use EEG as input
%
% Revision 1.96  2004/03/19 19:45:49  arno
% plotrad as number
%
% Revision 1.95  2004/03/19 19:45:16  arno
% conversion
%
% Revision 1.94  2004/03/19 19:04:53  scott
% help message and plotrad gui text
%
% Revision 1.93  2004/03/18 01:41:01  arno
% msgs
%
% Revision 1.92  2004/03/18 01:37:55  arno
% implemetning plotrad
%
% Revision 1.91  2004/03/18 00:28:41  arno
% debug shrink/skirt
%
% Revision 1.90  2004/03/08 19:24:49  arno
% shink -> shrink
%
% Revision 1.89  2004/02/13 19:13:56  arno
% nothing
%
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

function [chansout, chaninfo, urchans, com] = pop_chanedit(chans, orichaninfo, varargin);

urchans  = [];
com ='';
if nargin < 1
    help pop_chanedit;
    return;
end;
chansout = chans;
chaninfo   = [];
fig      = [];

if nargin < 2
    orichaninfo = [];
end;

if isempty(chans) | ~isnumeric(chans)
    % in case an EEG structure was given as input
    % -------------------------------------------
    if isfield(chans, 'chanlocs')
        dataset_input = 1;
        EEG           = chans;
        chans         = EEG(1).chanlocs;
        nchansori     = EEG.nbchan;
        if isfield(EEG, 'chaninfo')
            chaninfo = EEG(1).chaninfo;
        else chaninfo = [];
        end;
        if isfield(EEG, 'urchanlocs')
            urchans = EEG(1).urchanlocs;
        end;
    else
        nchansori     = 0;
        dataset_input = 0;
        chaninfo      = orichaninfo;
    end;

    % insert "no data channels" in channel structure
    % ----------------------------------------------
    nbchan = length(chans);
    [chans chaninfo] = insertchans(chans, chaninfo, nchansori);

    allfields = { 'labels' 'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' 'type' 'ref' 'urchan' 'datachan' };

    if isfield(chans, 'shrink')
        icadefs;
        if SHRINKWARNING
            warndlg2( [ 'You are currently shrinking channel locations for display.' 10 ...
                'A new option (more anatomically correct) is to plot channels' 10 ...
                'outside head limits so the shrink option has been disabled.' 10 ...
                '(Edit the icadefs file to disable this message)' ], 'Shrink factor warning');
        end;
        chans = rmfield(chans, 'shink');
    end;
    [chans shrinkorskirt plotrad] = checkchans(chans, allfields);

    % dealing with additional parameters
    % ----------------------------------
    if nargin > 1 & ~isstr(orichaninfo), % nothing
        if nargin > 2
            if ~isstr(varargin{1})
                urchans  = varargin{1};
                varargin = varargin(2:end);
            end;
        end;
    elseif nargin > 1
        varargin = { orichaninfo varargin{:} };
        clear chaninfo;
        chaninfo.shrink        = shrinkorskirt;
        chaninfo.plotrad       = plotrad;
    else
        chaninfo.shrink  = [];
        chaninfo.plotrad = [];
    end;

    nosevals       = { '+X' '-X' '+Y' '-Y' };
    if ~isfield(chaninfo, 'plotrad'), chaninfo.plotrad = []; end;
    if ~isfield(chaninfo, 'shrink'),  chaninfo.shrink = [];  end;
    if ~isfield(chaninfo, 'nosedir'), chaninfo.nosedir = nosevals{1}; end;
    oldchaninfo = chaninfo;
end;

if nargin < 3

    totaluserdat = {};
    % lookup channel locations if necessary
    % -------------------------------------
    if ~all(cellfun('isempty', {chans.labels})) && all(cellfun('isempty', {chans.theta}))
        [chans chaninfo urchans com] = pop_chanedit(chans, chaninfo, 'lookupgui', []);
        if ~isempty(com)
            totaluserdat = com;
            [chans chaninfo urchans com] = pop_chanedit(chans, chaninfo, com{:});
        end;
    end;

    commentfields = { 'Channel label ("label")', ...
                      'Polar angle ("theta")', 'Polar radius ("radius")', ...
                      'Cartesian X ("X")', ...
                      'Cartesian Y ("Y")', ...
                      'Cartesian Z ("Z")', ...
                      'Spherical horiz. angle ("sph_theta")', ...
                      'Spherical azimuth angle ("sph_phi")', ...
                      'Spherical radius ("sph_radius")' ...
                      'Channel type' 'Reference' ...
                      'Index in backup ''urchanlocs'' structure' ...
                      'Channel in data array (set=yes)' };

    % add field values
    % ----------------
    geometry = { 1 };
    tmpstr = sprintf('Channel information ("field_name"):');
    uilist = { { 'Style', 'text', 'string', tmpstr, 'fontweight', 'bold'  } };

    uiconvert = { ...
        { 'Style', 'pushbutton', 'string', 'Opt. head center', 'callback', 'pop_chanedit(gcbf, [], ''chancenter'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Rotate axis'     , 'callback', 'pop_chanedit(gcbf, [], ''forcelocs'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Transform axes'  , 'callback', 'pop_chanedit(gcbf, [], ''transform'', []);' } ...
        {  }, ...
        { 'Style', 'pushbutton', 'string', 'xyz -> polar & sph.', 'callback', 'pop_chanedit(gcbf, [], ''convert'', {''cart2all''});' }, ...
        { 'Style', 'pushbutton', 'string', 'sph. -> polar & xyz', 'callback', 'pop_chanedit(gcbf, [], ''convert'', {''sph2all'' });' }, ...
        { 'Style', 'pushbutton', 'string', 'polar -> sph. & xyz', 'callback', 'pop_chanedit(gcbf, [], ''convert'', {''topo2all''});' }, ...
        {  }, ...
        { 'Style', 'pushbutton', 'string', 'Set head radius', 'callback',   'pop_chanedit(gcbf, [], ''headrad'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Set channel types', 'callback', 'pop_chanedit(gcbf, [], ''settype'', []);' } ...
        { 'Style', 'pushbutton', 'string', 'Set reference', 'callback',     'pop_chanedit(gcbf, [], ''setref'' , []);' } ... 
        { } { } };

    % create text and edit for each field
    % -----------------------------------
    for index = 1:length(allfields)-1
        cbfield = [ 'valnumtmp   = str2num(get(findobj(gcbf, ''tag'', ''chaneditnumval''), ''string''));' ...
                    'pop_chanedit(gcbf, [], ''changefield'', { valnumtmp ''' allfields{index} ''' get(gcbo, ''string'') });' ...
                    'clear valnumtmp;' ];
        geometry = { geometry{:} [1.5 1 0.2 1] };
        uilist   = { uilist{:}, ...
            { 'Style', 'text', 'string', commentfields{index} }, ...
            { 'Style', 'edit', 'tag', [ 'chanedit' allfields{index} ], 'string', ...
            num2str(getfield(chans,{1}, allfields{index})), 'horizontalalignment', 'center', 'callback', cbfield } ...
              { } uiconvert{index} };
    end;

    % special checkbox for chandata field
    % -----------------------------------
    geometry = { geometry{:} [1.9 0.6 0.2 1] };
    cbfield = [ 'valnumtmp   = str2num(get(findobj(gcbf, ''tag'', ''chaneditnumval''), ''string''));' ...
                'pop_chanedit(gcbf, [], ''changefield'', { valnumtmp ''' allfields{end} ''' get(gcbo, ''value'') });' ...
                'clear valnumtmp;' ];
    uilist   = { uilist{:}, ...
        { 'Style', 'text', 'string', commentfields{end} }, ...
        { 'Style', 'checkbox', 'tag', [ 'chanedit' allfields{end}], 'string', '' 'value', 1 'callback', cbfield } { } uiconvert{end} };

    % add buttons
    % -----------
    geometry =  { geometry{:} [1] [1.15 0.5 0.6 1.9 0.4 0.4 1.15] [1.15 0.7 0.7 1 0.7 0.7 1.15] };
    cb_del    = [ 'valnum   = str2num(char(get(findobj(''tag'', ''chaneditnumval''), ''string'')));' ...
                  'pop_chanedit(gcbf, [], ''delete'', valnum);' ];
    cb_insert = [ 'valnum   = str2num(char(get(findobj(''tag'', ''chaneditnumval''), ''string'')));' ...
                  'pop_chanedit(gcbf, [], ''insert'', valnum);' ];
    cb_append = [ 'valnum   = str2num(char(get(findobj(''tag'', ''chaneditnumval''), ''string'')));' ...
                  'pop_chanedit(gcbf, [], ''append'', valnum);' ];

    uilist   = { uilist{:}, ...
        { }, ...
        { 'Style', 'pushbutton', 'string', 'Delete chan',  'callback', cb_del }, ...
        { },{ }, ...
        { 'Style', 'text'      , 'string', ['Channel number (of ' int2str(length(chans)) ')'], ...
        'fontweight', 'bold', 'tag', 'chaneditscantitle' }, { },{ },{ }, ...
        { 'Style', 'pushbutton', 'string', 'Insert chan',  'callback', cb_insert } ...
        { 'Style', 'pushbutton', 'string', '<<', 'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', -10);' ] } ...
        { 'Style', 'pushbutton', 'string', '<',  'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', -1);' ] } ...
        { 'Style', 'edit'      , 'string', '1', 'tag', 'chaneditnumval', 'callback', [ 'pop_chanedit(gcbf, []);' ] } ...
        { 'Style', 'pushbutton', 'string', '>',  'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', 1);' ] } ...
        { 'Style', 'pushbutton', 'string', '>>', 'callback', [ 'pop_chanedit(gcbf, [], ''movecursor'', 10);' ] } ...
        { 'Style', 'pushbutton', 'string', 'Append chan',  'callback', cb_append }, ...
        };

    % add sorting options
    % -------------------
    noseparam = strmatch(upper(chaninfo.nosedir), { '+X' '-X' '+Y' '-Y' });
    if isempty(noseparam), error('Wrong value for nose direction'); end;
    geometry = { geometry{:} [1] [0.9 1.3 0.6 1.1 0.9] [1] [1 1 1 1 1]};
    uilist   = { uilist{:},...
        { } ...
        { 'Style', 'pushbutton', 'string', 'Plot 2-D', 'callback', 'pop_chanedit(gcbf, [], ''plot2d'', []);' },...
        { 'Style', 'text', 'string', 'Plot radius (0.2-1, []=auto)'} ...
        { 'Style', 'edit', 'string', chaninfo.plotrad, 'tag', 'plotrad' 'callback' 'pop_chanedit(gcbf, [], ''plotrad'', []);' } ...
        { 'Style', 'popupmenu',  'string', 'Nose along +X|Nose along -X|Nose along +Y|Nose along -Y', ...
        'tag' 'nosedir' 'value',noseparam, 'callback' 'pop_chanedit(gcbf,[],''nosedir'',[]);' 'listboxtop' noseparam } ...
        { 'Style', 'pushbutton', 'string', 'Plot 3-D (xyz)',     'callback', 'pop_chanedit(gcbf, [], ''plot3d'', []);' } ...
        {}, ...
        { 'Style', 'pushbutton', 'string', 'Read locations',     'callback', 'pop_chanedit(gcbf,[],''load'',[]);' }, ...
        { 'Style', 'pushbutton', 'string', 'Read locs help',     'callback', 'pophelp(''readlocs.m'');' }, ...
        { 'Style', 'pushbutton', 'string', 'Look up locs',       'callback', 'pop_chanedit(gcbf,[], ''lookupgui'', []);' }, ...
        { 'Style', 'pushbutton', 'string', 'Save (as .ced)',     'callback', 'pop_chanedit(gcbf,[], ''save'',[]);' } ...
        { 'Style', 'pushbutton', 'string', 'Save (other types)'  'callback', 'pop_chanedit(gcbf,[], ''saveothers'',[]);' } ...
        };

    % evaluation of command below is required to center text (if
    % declared a text instead of edit, the uicontrol is not centered)
    comeval = [ 'set(findobj( ''tag'', ''chanediturchan''), ''style'', ''text'', ''backgroundcolor'', [.66 .76 1] );' ...
                'set(findobj( ''tag'', ''chaneditref''),    ''style'', ''text'', ''backgroundcolor'', [.66 .76 1] );' ...
                'set(findobj( ''tag'', ''ok''), ''callback'', ''pop_chanedit(gcbf, [], ''''return'''', []);'')' ];

    userdata.chans     = chans;
    userdata.nchansori = nchansori;
    userdata.chaninfo  = chaninfo;
    userdata.commands  = {};

    [results userdata returnmode] = inputgui( 'geometry', geometry, 'uilist', uilist, 'helpcom', ...
        'pophelp(''pop_chanedit'');', 'title', 'Edit channel info -- pop_chanedit()', ...
        'userdata', userdata, 'eval' , comeval );

    if length(results) == 0, 
        com = ''; 
        if dataset_input, chansout = EEG; end; 
        return; 
    end;

    % transfer events back from global workspace
    chans      = userdata.chans;
    chaninfo   = userdata.chaninfo;
    if ~isempty(userdata.commands)
        com = sprintf('%s=pop_chanedit(%s, %s);', inputname(1), inputname(1), vararg2str(userdata.commands));
    end;
else
    
    % call from command line or from a figure
    % ---------------------------------------
    currentpos = 0;
    if isnumeric(chans)
        fig         = chans;
        userdata    = get(fig, 'userdata');
        chans       = userdata.chans;
        nchansori   = userdata.nchansori;
        chaninfo    = userdata.chaninfo;
        currentpos  = str2num(get(findobj(fig, 'tag', 'chaneditnumval'), 'string'));
    end;
    
    args = varargin;
    % no interactive inputs
    % scan all the fields of g
    % ------------------------
    for curfield = 1:2:length(args)
        switch lower(args{curfield})
            case 'return'
                [tmpchans tmpfid] = getnodatchan(chans);
                if nchansori ~= 0 & nchansori ~= length(tmpchans)
                    if ~popask(strvcat(['The number of data channels (' int2str(length(tmpchans)) ') not including fiducials does not'], ...
                            ['correspond to the initial number of channels (' int2str(nchansori) '), so for consistency purposes'], ...
                            'new channel information will be ignored if this function was called from EEGLAB'))
                    else
                        set(findobj(fig, 'tag', 'ok'), 'userdata', 'stop');
                    end;
                else
                    set(findobj(fig, 'tag', 'ok'), 'userdata', 'stop');
                end;
                args = {};
            case 'plot3d', % GUI only
                tmpind = find(~cellfun('isempty', { chans.X }));
                if ~isempty(tmpind),
                    plotchans3d([ [ chans(tmpind).X ]' [ chans(tmpind).Y ]' [ chans(tmpind).Z ]'], { chans(tmpind).labels }); 
                else disp('cannot plot: no XYZ coordinates');
                end;
                args = {};
            case 'plot2d', % GUI only
                plotrad = str2num(get(findobj(fig, 'tag', 'plotrad'), 'string'));
                figure; topoplot([],chans, 'style', 'blank', 'drawaxis', 'on', 'electrodes', ...
                                 'labelpoint', 'plotrad', plotrad, 'chaninfo', chaninfo);
                args = {};
            case 'movecursor', % GUI only
                currentpos = max(1,min(currentpos+args{curfield+1},length(chans)));
                args = {};
            case 'plotrad',
                if isempty( args{curfield+1} )
                     args{curfield+1} = str2num(get(findobj(fig, 'tag', 'plotrad'), 'string'));
                end;
                chaninfo.plotrad = args{curfield+1};
            case 'forcelocs',
                if ~isempty(fig) % GUI BASED
                    [ comtmp tmpforce ] = forcelocs(chans); 
                    if ~isempty(tmpforce), 
                        args{curfield+1} = tmpforce{1};
                    end;
                end;
                if ~isempty(args{curfield+1})
                    chans = forcelocs(chans,args{curfield+1});
                    disp('Convert XYZ coordinates to spherical and polar');
                end;
            case 'chancenter',
                if ~isempty(fig)
                    [chans newcenter tmpcom] = pop_chancenter(chans);
                    args{curfield  } = 'eval';
                    args{curfield+1} = tmpcom;
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
                if ~isempty(fig) & ~strcmp(method, 'chancenter')
                    tmpButtonName=questdlg2( strvcat('This will modify fields in the channel structure', ...
                        'Are you sure you want to apply this function ?'), 'Confirmation', 'Cancel', 'Yes','Yes');
                    if ~strcmpi(tmpButtonName, 'Yes'), return; end;
                end;
                switch method
                    case 'chancenter',
                        if isempty(extraargs)
                            [X Y Z]=chancenter( [chans.X ]', [ chans.Y ]', [ chans.Z ]',[]);
                        else
                            [X Y Z]=chancenter( [chans.X ]', [ chans.Y ]', [ chans.Z ]', extraargs{:});
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
            case 'settype'                
                if ~isempty(fig)
                    args{curfield+1} = inputdlg2({'Channel indices' 'Type (e.g. EEG)' }, ...
                                             'Set channel type', 1, { '' '' }, 'pop_chanedit');
                end;
                try, tmpchans = args{curfield+1}{1}; tmptype = args{curfield+1}{2};catch, return; end;
                if isempty(tmpchans) & isempty(tmptype), return; end;
                if isstr(tmpchans)
                    tmpchans = eval( [ '[' tmpchans ']' ], 'settype: error in channel indices');
                end;
                if ~isstr(tmptype), tmptype = num2str(tmptype); end;
                for index = 1:length(tmpchans)
                    if tmpchans(index) > 0 & tmpchans(index) <= length(chans)
                        chans( tmpchans(index) ).type = tmptype;
                    end;
                end;
            case 'setref'
                if ~isempty(fig)
                    disp('Note that setting the reference only changes the reference labels');
                    disp('Use the re-referencing menu to change the reference');
                    args{curfield+1} = inputdlg2({'Channel indices' 'Reference (e.g. Cz)' }, ...
                                                 'Set channel reference', 1, { '' '' }, 'pop_chanedit');
                end;
                try, tmpchans = args{curfield+1}{1}; tmpref = args{curfield+1}{2};catch, return; end;
                if isempty(tmpchans) & isempty(tmpref), return; end;
                if isstr(tmpchans)
                    tmpchans = eval( [ '[' tmpchans ']' ], 'settype: error in channel indices');
                end;
                if ~isstr(tmpref), tmpref = num2str(tmpref); end;
                for index = 1:length(tmpchans)
                    if tmpchans(index) > 0 & tmpchans(index) <= length(chans)
                        chans( tmpchans(index) ).ref = tmpref;
                    end;
                end;
            case 'transform'
                if ~isempty(fig)
                    args{curfield+1} = inputdlg2({'Enter transform: (Ex: TMP=X; X=-Y; Y=TMP or Y(3) = X(2), etc.' }, ...
                                                  'Transform', 1, { '' }, 'pop_chanedit');
                end;
                try, tmpoper = args{curfield+1}; catch, return; end;
                if isempty(deblank(tmpoper)), return; end;
                if iscell(tmpoper), tmpoper = tmpoper{1}; end;
                tmpoper = [ tmpoper ';' ];
                [eloc, labels, theta, radius, indices] = readlocs(chans);
                if isempty(findstr(tmpoper, 'chans'))
                    try,
                        X          = [ chans(indices).X ];
                        Y          = [ chans(indices).Y ];
                        Z          = [ chans(indices).Z ];
                        sph_theta  = [ chans(indices).sph_theta ];
                        sph_phi    = [ chans(indices).sph_phi ];
                        sph_radius = [ chans(indices).sph_radius ];
                        eval(tmpoper);
                        
                        for ind = 1:length(indices)
                            chans(indices(ind)).X = X(min(length(X),ind));
                            chans(indices(ind)).Y = Y(min(length(Y),ind));
                            chans(indices(ind)).Z = Z(min(length(Z),ind));
                            chans(indices(ind)).theta  = theta(min(length(theta),ind));
                            chans(indices(ind)).radius = radius(min(length(radius),ind));
                            chans(indices(ind)).sph_theta  = sph_theta(min(length(sph_theta),ind));
                            chans(indices(ind)).sph_phi    = sph_phi(min(length(sph_phi),ind));
                            chans(indices(ind)).sph_radius = sph_radius(min(length(sph_radius),ind));
                        end;
                        
                        if     ~isempty(findstr(tmpoper, 'X')),          chans = convertlocs(chans, 'cart2all'); end;
                        if     ~isempty(findstr(tmpoper, 'Y')),          chans = convertlocs(chans, 'cart2all'); end;
                        if     ~isempty(findstr(tmpoper, 'Z')),          chans = convertlocs(chans, 'cart2all'); end;
                        if     ~isempty(findstr(tmpoper, 'sph_theta')),  chans = convertlocs(chans, 'sph2all');
                        elseif ~isempty(findstr(tmpoper, 'theta')),      chans = convertlocs(chans, 'topo2all'); end;
                        if     ~isempty(findstr(tmpoper, 'sph_phi')),    chans = convertlocs(chans, 'sph2all');  end;
                        if     ~isempty(findstr(tmpoper, 'sph_radius')), chans = convertlocs(chans, 'sph2all');
                        elseif ~isempty(findstr(tmpoper, 'radius')),     chans = convertlocs(chans, 'topo2all'); end;
                    catch, disp('Unknown error when applying transform'); end;
                else
                    eval(tmpoper);
                end;
                
            case 'headrad'
                if ~isempty(fig) % GUI
                     tmpres = inputdlg2({'Enter new head radius (same unit as DIPFIT head model):' }, ...
                                        'Head radius', 1, { '' }, 'pop_chanedit'); 
                     if ~isempty(tmpres),
                          args{ curfield+1 } = str2num(tmpres{1});
                     else return;
                     end;
                end;
                if ~isempty( args{ curfield+1 } )
                    allrad = [ chans.sph_radius ];
                    if length(unique(allrad)) == 1 % already spherical
                        chans = pop_chanedit(chans, 'transform', [ 'sph_radius = ' num2str( args{ curfield+1 } ) ';' ]);
                    else % non-spherical, finding best match
                        factor = args{ curfield+1 } / mean(allrad);
                        chans = pop_chanedit(chans, 'transform', [ 'sph_radius = sph_radius*' num2str( factor ) ';' ]);
                        disp('Warning: electrodes do not lie on a sphere. Sphere model fitting for');
                        disp('         dipole localization will work but generate many warnings');
                    end;
                    chans = convertlocs(chans, 'sph2all');
                end;
                
            case 'shrink'
                chans(1).shrink = args{ curfield+1 };
                
            case 'plotrad'
                chans(1).plotrad = args{ curfield+1 };
                
            case 'delete'
                chans(args{ curfield+1 })=[];
                
            case 'changefield'
                tmpargs = args{ curfield+1 };
                if length( tmpargs ) < 3
                    error('pop_chanedit: not enough arguments to change field value');
                end;
                eval([ 'chans(' int2str(tmpargs{1}) ').'  tmpargs{2} '=' reformat(tmpargs{3} ) ';' ]);
                
            case { 'insert' 'add' 'append' }
                tmpargs = args{ curfield+1 };
                allfields = fieldnames(chans);
                if isnumeric(tmpargs)
                    tmpargs2    = cell(1, length(allfields)+1);
                    tmpargs2{1} = tmpargs;
                    tmpargs     = tmpargs2;
                    if strcmpi(allfields{end}, 'datachan'), tmpargs{end} = 0; end;
                end;
                if length( tmpargs ) < length(allfields)+1
                    error('pop_chanedit: not enough arguments to change all field values');
                end;
                num = tmpargs{1};
                if strcmpi(lower(args{curfield}), 'append'), num=num+1; currentpos = currentpos+1; end;
                chans(end+1) = chans(end);
                chans(num+1:end) = chans(num:end-1);
                for index = 1:length( allfields )
                    chans = setfield(chans, {num}, allfields{index}, tmpargs{index+1});
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
                if ~isempty(fig) % GUI
                    [tmpf tmpp] = uigetfile('*.*', 'Load a channel location file'); 
                    drawnow;
                    if ~isequal(tmpf, 0),
                        tmpformats = readlocs('getinfos');
                        tmpformattype = { 'autodetect' tmpformats(1:end-1).type };
                        tmpformatstr  = { 'autodetect' tmpformats(1:end-1).typestring };
                        tmpformatdesc = { 'Autodetect file format from file extension' tmpformats(1:end-1).description };
                        %cb_listbox    = 'tmpdesc=get(gcbf, ''userdata''); set(findobj(gcbf, ''tag'', ''strdesc''), ''string'',  strmultiline([ ''File format: '' tmpdesc{get(gcbo, ''value'')} ], 30, 10)); clear tmpdesc;'' } }, ''pophelp(''''readlocs'''')'',' ...
                	    %        'Read electrode file'', tmpformatdesc, ''normal'', 4);
                        %txtgui = [ strmultiline([ 'File format: Autodetect file format from file extension'], 20, 10) 10 10 ];
                        %tmpfmt = inputgui( 'geometry', {[1 1]}, ...
                        %                   'uilist'  , { { 'style', 'text', 'string', txtgui 'tag' 'strdesc' }, ...
                        %                                 { 'style', 'listbox', 'string', strvcat(tmpformatstr) 'callback' '' } }, ...
                        %                   'geomvert', [10], ...
                        %                   'helpcom' , 'pophelp(''readlocs'');');
                        tmpfmt = inputgui( 'geometry', {[1 1 1] [1]}, ...
                                           'uilist'  , { { 'style', 'text', 'string', 'File format:' 'tag' 'strdesc' } {} {}, ...
                                                         { 'style', 'listbox', 'string', strvcat(tmpformatstr) 'callback' '' } }, ...
                                           'geomvert', [1 8], ...
                                           'helpcom' , 'pophelp(''readlocs'');');
                        if isempty(tmpfmt), 
                             args{ curfield+1 } = [];
                        else args{ curfield+1 } = { fullfile(tmpp, tmpf) 'filetype' tmpformattype{tmpfmt{1}} };
                        end;
                    else args{ curfield+1 } = [];
                	end;
                end;
                
                tmpargs = args{ curfield+1 };
                if ~isempty(tmpargs),
                    if isstr(tmpargs)
                        [chans] = readlocs(tmpargs);
                        chaninfo = [];
                        chaninfo.filename = tmpargs;
                    else
                        [chans] = readlocs(tmpargs{:});
                        chaninfo = [];
                        chaninfo.filename = tmpargs{1};
                    end;
                    for ind = 1:length(chans)
                        chans(ind).datachan = 1;
                    end;
                    
                    % backup file content etc...
                    % --------------------------
                    tmptext         = loadtxt( chaninfo.filename, 'delim', [], 'verbose', 'off', 'convert', 'off');
                    chaninfo.filecontent = strvcat(tmptext{:});
                    
                    % set urchan structure
                    % --------------------
                    urchans = chans;
                    for index = 1:length(chans)
                        chans(index).urchan = index;
                    end;
                end;
                
            case 'eval'
                tmpargs = args{ curfield+1 };
                eval(tmpargs);
                
            case 'saveothers'
                com = pop_writelocs(chans);
                args{ curfield }   = 'eval';
                args{ curfield+1 } = com;

            case 'save'
                if ~isempty(fig)
                     [tmpf tmpp] = uiputfile('*.ced', 'Save channel locs in EEGLAB .ced format');
                     drawnow;
                     args{ curfield+1 } = fullfile(tmpp, tmpf);
                end;
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
                
            case 'nosedir'
                nosevals = { '+X' '-X' '+Y' '-Y' };
                if ~isempty(fig)
                    tmpval = get(findobj(gcbf, 'tag', 'nosedir'), 'value');
                    args{ curfield+1 } = nosevals{tmpval};
                end;
                chaninfo.nosedir = args{ curfield+1 };
                if isempty(strmatch(chaninfo.nosedir, nosevals))
                    error('Wrong value for nose direction');
                end;
                
            case 'lookup'
                chaninfo.filename = args{ curfield+1 };
                if strcmpi(chaninfo.filename, 'standard-10-5-cap385.elp')
                    dipfitdefs;
                    chaninfo.filename = template_models(1).chanfile;
                elseif strcmpi(chaninfo.filename, 'standard_1005.elc')
                    dipfitdefs;
                    chaninfo.filename = template_models(2).chanfile;
                end;
                tmplocs = readlocs( chaninfo.filename, 'defaultelp', 'BESA' );
                [tmp ind1 ind2] = intersect(lower({ tmplocs.labels }), lower({ chans.labels }));
                if ~isempty(tmp)
                    chans = struct('labels', { chans.labels });
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
                        chans(ind2(index)).datachan   = 1;
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
                    if ~isfield(chans, 'type'), chans(1).type = []; end;
                end;
                if ~isempty(findstr(args{ curfield+1 }, 'standard_10')) & ...
                        ~isempty(findstr(args{ curfield+1 }, '.elc'))
                    chaninfo.nosedir = '+Y';
                else
                    chaninfo.nosedir = '+X';
                end;
                urchans = chans;
                for index = 1:length(chans)
                    chans(index).urchan = index;
                end;
                
            case 'lookupgui'
                standardchans = { 'Fp1' 'Fpz' 'Fp2' 'Nz' 'AF9' 'AF7' 'AF3' 'AFz' 'AF4' 'AF8' 'AF10' 'F9' 'F7' 'F5' ...
                    'F3' 'F1' 'Fz' 'F2' 'F4' 'F6' 'F8' 'F10' 'FT9' 'FT7' 'FC5' 'FC3' 'FC1' 'FCz' 'FC2' ...
                    'FC4' 'FC6' 'FT8' 'FT10' 'T9' 'T7' 'C5' 'C3' 'C1' 'Cz' 'C2' 'C4' 'C6' 'T8' 'T10' ...
                    'TP9' 'TP7' 'CP5' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4' 'CP6' 'TP8' 'TP10' 'P9' 'P7' 'P5' ...
                    'P3' 'P1' 'Pz' 'P2' 'P4' 'P6' 'P8' 'P10' 'PO9' 'PO7' 'PO3' 'POz' 'PO4' 'PO8' 'PO10' ...
                    'O1' 'Oz' 'O2' 'O9' 'O10' 'CB1' 'CB2' 'Iz' };
                [tmp1 ind1 ind2] = intersect( lower(standardchans), lower({ chans.labels }));
                if ~isempty(tmp1) | isfield(chans, 'theta')
                    
                    % finding template location files
                    % -------------------------------
                    setmodel = [ 'tmpdat = get(gcbf, ''userdata'');' ...
                        'tmpval = get(gcbo, ''value'');' ...
                        'set(findobj(gcbf, ''tag'', ''elec''), ''string'', tmpdat{tmpval});' ...
                        'clear tmpval tmpdat;' ];
                    try
                        EEG = eeg_emptyset; % for dipfitdefs
                        dipfitdefs;
                        userdatatmp = { template_models(1).chanfile template_models(2).chanfile };
                        clear EEG;
                    catch, userdatatmp = { 'Standard-10-5-Cap385.sfp' 'Standard-10-5-Cap385.sfp' };
                    end;
                    
                    % other commands for help/load
                    % ----------------------------
                    comhelp = [ 'warndlg2(strvcat(''The template file depends on the model'',' ...
                        '''you intend to use for dipole fitting. The default file is fine for'',' ...
                        '''spherical model.'');' ];
                    commandload = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
                        'if filename ~=0,' ...
                        '   set(findobj(''parent'', gcbf, ''tag'', ''elec''), ''string'', [ filepath filename ]);' ...
                        'end;' ...
                        'clear filename filepath tagtest;' ];
                    if ~isfield(chans, 'theta'),                    message =1;
                    elseif all(cellfun('isempty', {chans.theta })), message =1;
                    else                                            message =2;
                    end;
                    if message == 1
                        textcomment = strvcat('Only channel labels are present currently, but some of these labels have known', ...
                            'positions. Do you want to look up coordinates for these channels using the electrode', ...
                            'file below? If you have a channel location file for this dataset, press cancel, then', ...
                            'use button "Read location" in the following gui. If you do not know, just press OK.');
                    else
                        textcomment = strvcat('Some channel labels may have known locations.', ...
                            'Do you want to look up coordinates for these channels using the electrode', ...
                            'file below? If you do not know, press OK.');
                    end;
                    uilist = { { 'style' 'text' 'string' textcomment } ...
                        { 'style' 'popupmenu'  'string' [ 'use BESA file for 4-shell dipfit spherical model' ...
                        '|use MNI coordinate file for BEM dipfit model' ] ...
                        'callback' setmodel } ...
                        { } ...
                        { 'style' 'edit'       'string' userdatatmp{1} 'tag' 'elec' } ...
                        { 'style' 'pushbutton' 'string' '...' 'callback' commandload } };
                    
                    res = inputgui( { 1 [1 0.3] [1 0.3] }, uilist, 'pophelp(''pop_chanedit'')', 'Look up channel locations?', userdatatmp, 'normal', [4 1 1] );
                    if ~isempty(res)
                        [chans chaninfo] = pop_chanedit(chans, chaninfo, 'lookup', res{2} );
                        args = { 'lookup', res{2} };
                    else
                        args = {};
                    end;
                end;
        end;
    end;
    
end;

% call from a figure
% ------------------
if ~isempty(fig)
    userdata.chans    = chans;
    userdata.chaninfo = chaninfo;
    userdata.commands = { userdata.commands{:} args{:} };
    set(fig, 'userdata', userdata);

    set(findobj(fig, 'tag', 'chaneditnumval'), 'string', num2str(currentpos));
    set(findobj(fig, 'tag', 'chaneditscantitle'), 'string', ['Channel number (of ' int2str(length(chans)) ')']);

    % update GUI with current channel info
    allfields = fieldnames(chans);
    for index = 1:length(allfields)
        obj = findobj(fig, 'tag', [ 'chanedit' allfields{index}]);
        if strcmpi(allfields{index}, 'datachan')
            set(obj, 'value', getfield(chans(currentpos), allfields{index}));
        else
            set(obj, 'string', num2str(getfield(chans(currentpos), allfields{index})));
        end;
    end;
else
    % command line call or OK in figure
    % ---------------------------------
    if isfield(chans, 'sph_phi_besa'  ), chans = rmfield(chans, 'sph_phi_besa'); end;
    if isfield(chans, 'sph_theta_besa'), chans = rmfield(chans, 'sph_theta_besa'); end;

    % move no data channels to info structure
    % ---------------------------------------
    [chans chaninfo.nodatchans] = getnodatchan(chans);
    if isempty(chaninfo.nodatchans), chaninfo = rmfield(chaninfo, 'nodatchans'); end;

    if dataset_input,
         if nchansori == length(chans)
             for index = 1:length(EEG)
                 EEG(index).chanlocs = chans;
                 EEG(index).chaninfo = chaninfo;
             end;
         else
             disp('Wrong channel structure size, changes ignored');
         end;
         chansout = EEG;
    else chansout = chans;
    end;
end;

return;

% format the output field
% -----------------------
function strval = reformat( val )
if isnumeric(val) & isempty(val), val = '[]'; end;
if isstr(val), strval = [ '''' val '''' ];
else           strval = num2str(val);
end;

% extract text using tokens (not used)
% ------------------------------------
function txt = inserttxt( txt, tokins, tokfind);
locfind = findstr(txt, tokfind);
for index = length(locfind):-1:1
    txt = [txt(1:locfind(index)-1) tokins txt(locfind(index):end)];
end;

% ask for confirmation
% --------------------
function num = popask( text )
ButtonName=questdlg2( text, ...
    'Confirmation', 'Cancel', 'Yes','Yes');
switch lower(ButtonName),
    case 'cancel', num = 0;
    case 'yes',    num = 1;
end;

% shrink and skirt factors
% ------------------------
function [chans, shrinkorskirt, plotrad]= checkchans(chans, fields);

shrinkorskirt = [];
plotrad  = [];
if isfield(chans, 'plotrad'),
    plotrad = chans(1).plotrad;
    chans = rmfield(chans, 'plotrad');
    if isstr(plotrad) & ~isempty(str2num(plotrad)), plotrad = str2num(plotrad); end;
end;
if isfield(chans, 'shrink')
    shrinkorskirt = 1;
    if ~isstr(chans(1).shrink)
        plotrad = 0.5/(1-chans(1).shrink); % convert old values
    end;
    chans = rmfield(chans, 'shrink');
end;

for index = 1:length(fields)
    if ~isfield(chans, fields{index})
        if ~strcmpi(fields{index}, 'datachan')
            chans = setfield(chans, {1}, fields{index}, []);
        elseif ~strcmpi(fields{index}, 'type')
            for indchan = 1:length(chans)
                chans = setfield(chans, {indchan}, num2str(fields{index}), indchan);
            end;
        else
            for indchan = 1:length(chans)
                chans = setfield(chans, {indchan}, fields{index}, indchan);
            end;
        end;
    end;
end;
if exist('orderfields') == 2
    try,
        chans = orderfields(chans, fields);
    catch, end;
end;

% separate data channels from non-data channels
% ---------------------------------------------
function [chans, fids] = getnodatchan(chans)
if isfield(chans, 'datachan')
    for ind = 1:length(chans)
        if isempty(chans(ind).datachan)
            chans(ind).datachan = 0;
        end;
    end;
    alldatchans = [ chans.datachan ];
    fids  = chans(find(alldatchans == 0));
    chans = chans(find(alldatchans));
else
    fids = [];
end;

% fuse data channels and non-data channels
% ----------------------------------------
function [chans, chaninfo] = insertchans(chans, chaninfo, nchans)
if nargin < 3, nchans = length(chans); end;
for ind = 1:length(chans)
    chans(ind).datachan = 1;
end;
if length(chans) > nchans % reference at the end of the structure
    chans(end).datachan = 0;
end;
if isfield(chaninfo, 'nodatchans')
    if ~isempty(chaninfo.nodatchans)
        chanlen = length(chans);
        for index = 1:length(chaninfo.nodatchans)
            fields = fieldnames( chaninfo.nodatchans );
            ind = chanlen+index;
            for f = 1:length( fields )
                chans = setfield(chans, { ind }, fields{f}, getfield( chaninfo.nodatchans, { index },  fields{f}));
            end;
            chans(ind).datachan = 0;
        end;
        disp('Non-data channels have been added to the end of the channel structure');
        chaninfo = rmfield(chaninfo, 'nodatchans');
        
        % put these channels first
        % ------------------------
        % tmp = chans(chanlen+1:end);
        % chans(length(tmp)+1:end) = chans(1:end-length(tmp));
        % chans(1:length(tmp)) = tmp;
    end;
end;

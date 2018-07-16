% pop_topochansel() - pop up a topographic interface to select channels
%
% Pops up a topographic interface to select multiple channels with the
% mouse. Click a polygon around the electrodes you wish to select. Right
% click to finish.
% 
% Usage:
%   >> [chanlist] = pop_topochansel(chanlocs,selection);
%
% Inputs:
%   chanstruct     - channel structure. See readlocs()
%                    (optional if EEG structure with chanlocs in caller or
%                    base workspace)
%   selection      - currently selected electrodes.
%                    (optional)
%
%
% Output:
%   chanlist      - indices of selected channels
%   cellchannames - names of selected channel names in a cell array
%   strchannames  - names of selected channel names in a concatenated string
%                   (channel names are separated by space characters)
%

function [chanlist, cellchannames, strchannames] = pop_topochansel(chanlocs,select,varargin)
if nargin == 0 || isempty(chanlocs)
    try
        chanlocs = evalin('caller','chanlocs');
    catch
        try
            tmp = evalin('caller','EEG');
            chanlocs = tmp.chanlocs;
        catch
            error('tell me where the electrodes are'),
        end
    end
end
if nargin >= 2
    if ischar(select)|| iscellstr(select)
        select = chnb(select,{chanlocs.labels});
    end
end
g = finputcheck( varargin, ...
                 { 'labels'    'string'    {'on' 'off'}   'off'
                  'cellstrout' 'string'    {'on' 'off'}   'off'
                    },'ignore');

        
h0 = figure(23537);clf;
set(gcf,'numbertitle','off','name','Channel selection');
topoplot([],chanlocs,'electrodes','on', 'style','blank');

chi = get(gca,'children');
% chi = chi(1:numel(chanlocs));
todel = [];
for i = 1:numel(chi)
    try
        if numel(get(chi(i),'XData')) == sum(~emptycells({chanlocs.X}))
            continue % we found electrodes
        else
            todel(end+1) = i;
        end
    catch
        todel(end+1) = i;
    end
end
chi(todel) = [];
hold on
X = get(chi,'XData');
Y = get(chi,'YData');
delete(chi)
if strcmp(g.labels,'on')
    for i = 1:numel(X)
        text(X(i),Y(i),chanlocs(i).labels);
    end
end

if exist('select','var')
    plot(X(select),Y(select),'r.','markersize',20);
else
    select = [];
end
[dum dum chanlist] = lasso(X,Y);
if isempty(chanlist)
    disp('No new elecs selected. Keeping old selection')
    chanlist = select;
end
strchannames = '';
for i = 1:numel(chanlist)
    strchannames = [strchannames chanlocs(chanlist(i)).labels];
    if i ~= numel(chanlist)
        strchannames = [strchannames ' '];
    end
    cellchannames{i} = chanlocs(chanlist(i)).labels;
end
if strcmp(g.cellstrout,'on')
    chanlist = cellchannames;
end
if nargout == 0
    disp(strchannames)
    disp(chanlist)
    clear
end
% close(h0);


function [selx,sely,indexnr]=lasso(x,y)

% lasso -  enables the selection/encircling of (clusters of) events in a scatter plot by hand
%          using the mouse
%
% Input:    x,y                 - a set of points in 2 column vectors.
% Output:   selx,sely,indexnr   - a set of selected points in 3 column vectors
%
% Note:   After the scatter plot is given, selection by mouse is started after any key press.
%         This is done to be able to ZOOM or CHANGE AXES etc. in the representation before selection
%         by mouse.
%         Encircling is done by pressing subsequently the LEFT button mouse at the requested positions
%         in a scatter plot.
%         Closing the loop is done by a RIGHT button press.
%
% T.Rutten V2.0/9/2003
% downloaded and adapted from mathworks fileexchange by M. Chaumon 2011

plot(x,y,'ob')

las_x=[];
las_y=[];

c=1;

key=0;

while c==1
    try
        [a,b,c]=ginput(1);
        las_x=[las_x;a];las_y=[las_y;b];
        line(las_x,las_y)
    catch
        selx = [];sely = []; indexnr = [];
        return
    end
end

las_x(length(las_x)+1)=las_x(1);
las_y(length(las_y)+1)=las_y(1);

line(las_x,las_y)
pause(.2)

in=inpolygon(x,y,las_x,las_y);

ev_in=find(in>0);

selx=x(ev_in);
sely=y(ev_in);
plot(selx,sely,'r.','markersize',20);
drawnow

indexnr=ev_in;

function nb = chnb(channame, varargin)

% chnb() - return channel number corresponding to channel names in an EEG
%           structure
%
% Usage:
%   >> [nb] = chnb(channame);
%   >> [nb] = chnb(channame, labels);
%
% Input:
%   channame      - name of channels to search. Either a string with space
%                   separated channel names, or a cell array of strings.
%                   Note that regular expressions can be used to match
%                   several channels. See regexp.
%   labels        - channel labels as found in EEG.chanlocs.labels.
%
% Output:
%   nb            - channel numbers in the EEG structure found in the
%                   caller workspace (i.e. where the function is called
%                   from) or in the base workspace, if no EEG structure
%                   exists in the caller workspace.
%
error(nargchk(1,2,nargin));
if nargin == 2
    labels = varargin{1};
else
    
    try
        EEG = evalin('caller','EEG');
    catch
        try
            EEG = evalin('base','EEG');
        catch
            error('Could not find EEG structure');
        end
    end
    if not(isfield(EEG,'chanlocs'))
        error('No channel list found');
    end
    labels = {EEG.chanlocs.labels};
end
if ischar(channame)
    tmp = regexp(channame,'(\S*) ?','tokens');
    channame = {};
    for i = 1:numel(tmp)
        channame{i} = tmp{i}{1};
    end
    if isempty(channame)
        nb = [];
        return
    end
end

nb = regexpcell(labels,channame,'exactignorecase');


function idx = regexpcell(c,pat, cmds)

% idx = regexpcell(c,pat, cmds)
%
% Return indices idx of cells in c that match pattern(s) pat (regular expression).
% Pattern pat can be char or cellstr. In the later case regexpcell returns
% indexes of cells that match any pattern in pat.
%
% cmds is a string that can contain one or several of these commands:
% 'inv' return indexes that do not match the pattern.
% 'ignorecase' will use regexpi instead of regexp
% 'exact' performs an exact match (regular expression should match the whole strings in c).
% 'all' (default) returns all indices, including repeats (if several pat match a single cell in c).
% 'unique' will return unique sorted indices.
% 'intersect' will return only indices in c that match ALL the patterns in pat.
% 
% v1 Maximilien Chaumon 01/05/09
% v1.1 Maximilien Chaumon 24/05/09 - added ignorecase
% v2 Maximilien Chaumon 02/03/2010 changed input method.
%       inv,ignorecase,exact,combine are replaced by cmds

error(nargchk(2,3,nargin))
if not(iscellstr(c))
    error('input c must be a cell array of strings');
end
if nargin == 2
    cmds = '';
end
if not(isempty(regexpi(cmds,'inv', 'once' )))
    inv = true;
else
    inv = false;
end
if not(isempty(regexpi(cmds,'ignorecase', 'once' )))
    ignorecase = true;
else
    ignorecase = false;
end
if not(isempty(regexpi(cmds,'exact', 'once' )))
    exact = true;
else
    exact = false;
end
if not(isempty(regexpi(cmds,'unique', 'once' )))
    combine = 2;
elseif not(isempty(regexpi(cmds,'intersect', 'once' )))
    combine = 3;
else
    combine = 1;
end

if ischar(pat)
    pat = cellstr(pat);
end

if exact
    for i_pat = 1:numel(pat)
        pat{i_pat} = ['^' pat{i_pat} '$'];
    end
end
    
for i_pat = 1:length(pat)
    if ignorecase
        trouv = regexpi(c,pat{i_pat}); % apply regexp on each pattern
    else
        trouv = regexp(c,pat{i_pat}); % apply regexp on each pattern
    end
    idx{i_pat} = [];
    for i = 1:numel(trouv)
        if not(isempty(trouv{i}))% if there is a match, store index
            idx{i_pat}(end+1) = i;
        end
    end
end
switch combine
    case 1
        idx = [idx{:}];
    case 2
        idx = unique([idx{:}]);
    case 3
        for i_pat = 2:length(pat)
            idx{1} = intersect(idx{1},idx{i_pat});
        end
        idx = idx{1};
end
if inv % if we want to invert result, then do so.
    others = 1:numel(trouv);
    others(idx) = [];
    idx = others;
end

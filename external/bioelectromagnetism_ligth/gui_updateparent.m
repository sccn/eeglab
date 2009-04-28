function [p] = gui_updateparent(UserData,focus)

% gui_updateparent - General GUI data handing for EEG Toolbox
%
% Usage: [p] = gui_updateparent(UserData,focus)
%
% 'UserData' is obtained from get(gcbf,'UserData'), where all
% the necessary parameters and stored for each GUI, including
% the handles of any GUI parent.
%
% 'focus' is a boolean option that switches whether the 
% parent is the current figure (1, default) or not (0).
% 
% This function returns the p struct to both the matlab
% workspace and any direct parent of the GUI calling 
% this function.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:56 $

% Licence:  GNU GPL, no express or implied warranties
% History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
%                    - extracted from individual gui* functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('focus','var'), focus = 1; end

if isfield(UserData,'parent'),
    if isfield(UserData.parent,'gui'),
        % Get the userdata from the parent
        parent = get(UserData.parent.gui,'UserData');
        if isfield(parent,'p') & isfield(UserData,'p'),
            % Update the parent p structure
            parent.p = UserData.p;
            set(UserData.parent.gui,'UserData',parent);
            % Make the parent the current figure
        end
        if focus,
            figure(UserData.parent.gui); end
    end
end

if isfield(UserData,'p'),
    if ~isempty(UserData.p),
       [p] = UserData.p;
    end
end

return

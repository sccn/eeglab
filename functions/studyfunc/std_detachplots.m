% std_detachplots() - Given a figure with subplots and several lines per axis will add a callback to
% each axis specified in the 'figtitles' input. The callback consist in a
% figure with all the detached individuals lines.
%
% Usage:
%   >>   std_detachplots('','','data',data 'figtitles', alltitlestmp,'sbtitles',sbtitles,'handles', handles);
%
% Inputs:
%      data        - 
%      g.timevec   - vector of loaded EEG datasets
%
% Optional inputs:
%
% Outputs:
%
% See also:
%   std_plotinfocluster
%
% Author: Ramon Martinez-Cancino, SCCN, 2014
%
% Copyright (C) 2014  Ramon Martinez-Cancino,INC, SCCN
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
function std_detachplots(hObject,eventdata,varargin)

% display help if not enough arguments
if nargin < 2
    help std_detachplots;
    return;
end
icadefs;
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end;
catch
    disp('std_detachplots() error: calling convention {''key'', value, ... } error'); return;
end;

try g.data;         catch, g.data          =   [];  end; % Name of plots
try g.figtitles;    catch, g.figtitles     =   [];  end; % Name of plots
try g.sbtitles;     catch, g.sbtitles      =   [];  end; % Name of plots
try g.handles;      catch, g.handles       =   [];  end; % Handles of figure
try g.flagstd;      catch, g.flagstd       =   1;   end; % Plot std band around mean
try g.xlabel;       catch, g.xlabel        =   '';  end; % xlabel
try g.ylabel;       catch, g.ylabel        =   '';  end; % ylabel
try g.timevec;      catch, g.timevec       =   '';  end; % Time or freq vector

% Checking data
 if isempty(g.handles) && any([isempty(g.data) isempty(g.figtitles)])
     error('std_detachplots : Check entries. Options ''handles'', ''data'' and ''figtitles'' must be provided');
 end

if iscell(g.data)
    nplots = numel(g.data(:));
else
    nplots = 1 ;
    g.data = {g.data};
end

%  Checking sbtitles 
if isempty(g.sbtitles)
    for i = 1:numel(g.data)
        c = 1;
        if ~isempty(g.data{i})
            for j = 1:size(g.data{i},2)
                g.sbtitles{i}{j} = {['Line ' num2str(c)]};
                c = c+1;
            end
%         else
%             g.sbtitles{i} = [];
        end
    end
end

% Plot goes here
%--------------------------------------------------------------------------
if isempty(g.handles) 
    for i_nplots = 1 : nplots
        idata = g.data{i_nplots};
        len = size(idata,2);
        if len > 0 % A non-empty cluster
            % Getting the mean
            meandata = mean(idata,2);
            stddata  = std(idata,0,2);
            lower    = meandata-stddata;
            upper    = meandata+stddata;
            
            hplot = figure('name', g.figtitles, 'NumberTitle','off');
            rowcols(2) = ceil(sqrt(len + 4));
            rowcols(1) = ceil((len+4)/rowcols(2));
            
            for k = 1:len 
                %--- first sbplot row ----
                if k <= rowcols(2) - 2
                    figure(hplot);
                    sbplot(rowcols(1),rowcols(2),k+2);
                    hold on;
                    plotlines(k,idata,meandata,lower, upper,g.sbtitles{k},g);
                else 
                    figure(hplot)
                    sbplot(rowcols(1),rowcols(2),k+4);
                    hold on;
                    plotlines(k,idata,meandata,lower, upper,g.sbtitles{k},g);
                end
            end
            % Plot all figure
            figure(hplot)
            sbplot(rowcols(1),rowcols(2),[1 rowcols(2)+2 ]);
            hold on;
            plotlines(1:len,idata,meandata,lower,upper,g.figtitles,g);
            set(gcf,'Color', BACKCOLOR);
            orient tall  % fill the figure page for printing
        end % Finished one cluster plot
    end
else
    xlabelval = '';
    ylabelval = '';
    
    % Match Children handles based on titles provided 
    for i = 1: nplots
        htemp          = findall(g.handles,'String', g.figtitles{i});
        if all([~isempty(htemp) ~isempty(g.data(i))])
        handlestemp{i} = htemp(1);  
        % Getting x label from handles
        if isempty(xlabelval)
            xlabelval = handlestemp{i}.Parent.XLabel.String;
        end
        % Getting y label from handles
        if isempty(ylabelval)
            ylabelval = handlestemp{i}.Parent.YLabel.String;
        end
        % Getting timevec from handles
        if isempty(g.timevec)
            g.timevec = handlestemp{i}.Parent.Children(i).XData;
        end
        else
            handlestemp{i} = [];
        end
    end
    c = 0;
    for i = 1: nplots
        if ~isempty(handlestemp{i})
            c = c + 1;
            % For Axis
            set( handlestemp{i}.Parent, 'ButtonDownFcn',{@std_detachplots,'data',g.data{i},'timevec',g.timevec,'handles' ,'',...
                'figtitles',g.figtitles{i},...
                'sbtitles',g.sbtitles{c},...
                'xlabel',xlabelval,...
                'ylabel',ylabelval,...
                });
            % For lines
            set( handlestemp{i}.Parent.Children, 'ButtonDownFcn',{@std_detachplots,'data',g.data{i},'timevec',g.timevec,'handles'  ,'',...
                'figtitles',g.figtitles{i},...
                'sbtitles',g.sbtitles{c},...
                'xlabel',xlabelval,...
                'ylabel',ylabelval,...
                });
        end
    end
end

function plotlines(kindx,idata,meandata,lower,upper,sbtitle,g)
if g.flagstd                           % plot std band
    eeglabciplot(lower,upper,g.timevec, 'r', 0.2);
    axis tight;
end
for i = 1:length(kindx)
    plot(g.timevec,idata(:,kindx(i)),'b','LineWidth', 0.1); % plot g.data
end
plot(g.timevec,meandata,'r','LineWidth', 0.1);              % plot mean
if length(kindx)>1
    xlabel(g.xlabel);
    ylabel(g.ylabel);
end
box on;
grid on;
axis tight;
title(sbtitle, 'interpreter', 'none');

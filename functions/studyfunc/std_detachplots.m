% std_detachplots() -  Given a figure with subplots and several lines per axis, will add a callback to
% each axis specified in the 'figtitles' input. The callback consist in a figure with all the detached
% individuals lines.
%
% Usage:
%   >>   std_detachplots('','','data',data 'figtitles', alltitlestmp,'sbtitles',sbtitles,'handles', handles);
%
% Inputs:
%      data        - Cell array containing the data matrices for each plot in the same order showed in the figure
%      figtitles   - Cell array of the titles of each individual axes in
%                    the  figure. The titles must correspond. The function
%                    use this value to find the right hanlde of the axis
%      sbtitles    - Cell array of cell arrays with the titles for each
%                    detached line per axis. i.e. {{'Axis1 line1' 'Axis1 line2'} {'Axis2 line1' 'Axis2 line2'}}
%      handles     - Handles of the main figure who contain all the
%                    subplots
%      flagstd     - Flag to plot the Standar Deviation  {default: 1} means 'on'
%
% See also:
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
try g.filter;       catch, g.filter        =   '';  end; % Low pass filter freq

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
        end
    end
end

% Plot goes here
%--------------------------------------------------------------------------
if isempty(g.handles)
    for i_nplots = 1 : nplots
        idata = g.data{i_nplots};
        
        % Filtering data to be plotted
        if ~isempty(g.filter), idata = myfilt(idata, 1000/(g.timevec(2)-g.timevec(1)), 0, g.filter); end;
        
        len = size(idata,2);
        if len > 0 % A non-empty cluster
            % Getting the mean
            meandata = mean(idata,2);
            stddata  = std(idata,0,2);
            
            if license('checkout', 'statistics_toolbox')
                SEM = std(idata,0,2)/sqrt(size(idata,2));       % Standard Error
                ts = tinv([0.025  0.975],size(idata,2)-1);      % T-Score
                lower    = meandata + ts(1)*SEM;                % CI
                upper    = meandata + ts(2)*SEM;                % CI
            else
                lower    = meandata-2*stddata;
                upper    = meandata+2*stddata;
            end
                
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
            orient tall;
        end
    end
else
    xlabelval = '';
    ylabelval = '';
    % Match Children handles based on titles provided
    c = 0;
    for i = 1: nplots
        htemp          = findall(g.handles,'String', g.figtitles{i});
        if all([~isempty(htemp) ~isempty(g.data(i))])
            handlestemp{i} = htemp(1);
            % Getting x label from handles
            if isempty(xlabelval)
                xlabelval = get(get(get(handlestemp{i},'Parent'),'Xlabel'),'String'); % xlabelval = handlestemp{i}.Parent.XLabel.String;
            end
            % Getting y label from handles
            if isempty(ylabelval)
                ylabelval = get(get(get(handlestemp{i},'Parent'),'Ylabel'),'String'); % ylabelval = handlestemp{i}.Parent.YLabel.String;
            end
            % Getting timevec from handles
            if isempty(g.timevec)
                tmp = get(get(get(handlestemp{i},'Parent'),'Children'));
                g.timevec = tmp(i).XData; % g.timevec = handlestemp{i}.Parent.Children(i).XData;
            end
        else
            handlestemp{i} = [];
        end
        
        % Callback setting
        if ~isempty(handlestemp{i})
            c = c + 1;
            % For Axis
            set( get( handlestemp{i}, 'Parent'), 'ButtonDownFcn',{@std_detachplots,...
                'data'     ,g.data{i},...
                'timevec'  ,g.timevec,...
                'handles'  ,'',...
                'figtitles',g.figtitles{i},...
                'sbtitles' ,g.sbtitles{c},...
                'xlabel'   ,xlabelval,...
                'ylabel'   ,ylabelval,...
                'filter'   , g.filter,...
                });
            % For lines
            set( get(get( handlestemp{i}, 'Parent'), 'Children'), 'ButtonDownFcn',{@std_detachplots,...
                'data'     ,g.data{i},...
                'timevec'  ,g.timevec,...
                'handles'  ,'',...
                'figtitles',g.figtitles{i},...
                'sbtitles' ,g.sbtitles{c},...
                'xlabel'   ,xlabelval,...
                'ylabel'   ,ylabelval,...
                'filter'   , g.filter,...
                });
        end
    end
end

function plotlines(kindx,idata,meandata,lower,upper,sbtitle,g)
for i = 1:length(kindx)
    plot(g.timevec,idata(:,kindx(i)),'b','LineWidth', 0.1); % plot g.data
end
plot(g.timevec,meandata,'r','LineWidth', 0.1);              % plot mean

if g.flagstd                                                % plot std band
    eeglabciplot(lower,upper,g.timevec, 'r', 0.2);
    axis tight;
end

if length(kindx)>1
    xlabel(g.xlabel);
    ylabel(g.ylabel);
end
box on;
grid on;
axis tight;
title(sbtitle, 'interpreter', 'none');

% rapid filtering for ERP (from std_plotcurve)
% -----------------------
function tmpdata2 = myfilt(tmpdata, srate, lowpass, highpass); 
    bscorrect = 1;
    if bscorrect
        % Getting initial baseline
        bs_val1  =  mean(tmpdata,1);
        bs1      = repmat(bs_val1, size(tmpdata,1), 1);
    end
    
    % Filtering
    tmpdata2 = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3)*size(tmpdata,4));
    tmpdata2 = eegfiltfft(tmpdata2',srate, lowpass, highpass)';
    tmpdata2 = reshape(tmpdata2, size(tmpdata,1), size(tmpdata,2), size(tmpdata,3), size(tmpdata,4));
    
    if bscorrect
        % Getting after-filter baseline
        bs_val2  =  mean(tmpdata2,1);
        bs2      = repmat(bs_val2, size(tmpdata2,1), 1);
        
        % Correcting the baseline
        realbs = bs1-bs2;
        tmpdata2 = tmpdata2 + realbs;
    end
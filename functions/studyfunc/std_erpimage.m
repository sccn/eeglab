% subjerpimage() - Plot an erpimage which is a mixture of erpimage
%                  for each subject
%
% Usage:
%   >> subjerpimage( subjind, data, sortvar, times, titleim, movewin, ...
%                    decim, erpimage_args ...);
%   >> subjerpimage( subjind, data, sortvar, times, titleim, opterp1, opterp2);
%
% Inputs:
%  subjind       - index of each subject for all trials
%  data          - data
%  sortvar       - sorting variable
%  times         - time vector (only forwarded to erpimage())
%  titleim       - image title (only forwarded to erpimage())
%  movewin       - moving average. Can have 2 values, the first one for
%                  inter-subject smoothing , the second one for across 
%                  subject smooting. For instance [10 2] will apply a
%                  10-trial smoothing average across subjects and 2-trial
%                  smoothing average across subjects. A negative value
%                  will subtract the smoothing average from the erpimage
%                  with a smoothing average 4-times smaller.
%  decim         - decimation factor. Enter the number of trials you want
%                  for each subject. Can have 2 values, the first one for
%                  inter-subject decimation , the second one for across 
%                  subject decimation.
%  erpimage_args - other erpimage() arguments. 
%  opterp1       - cell array of erp image argument for within subject 
%                  erpimage()
%  opterp2       - cell array of erp image argument for accross subject 
%                  erpimage()
%
% Outputs:
%  allphases   - all phase from erpimage
%  allerpimage - all sorting variables from erpimage
%  subjamp     - value of amplitude contribution [0 to 1] for each trial.
%                1 indicates that all subject contribute equally to the
%                amplitude of each trial. 0 indicates that only one
%                subject dominates all trials' amplitude.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2004

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Arnaud Delorme
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

%$Log: not supported by cvs2svn $
%Revision 1.2  2005/02/15 02:21:43  arno
%plot_indiv option
%
%Revision 1.1  2005/02/03 00:10:48  arno
%Initial revision
%

function [allphases, allsortvar, subjamptime, subjamptrial, globalent ] = subjerpimage( subjind, data, sortvar, times, titleim, movewin, decim, varargin );
    
    plot_indiv = 0; % 0 or 1
    
    if nargin < 7
        help subjerpimage
        return;
    end;
    if iscell(movewin)
        erparg1 = movewin;
        erparg2 = decim;
        movewin = [ erparg1{1} erparg2{1} ];
        decim   = [ erparg1{2} erparg2{2} ];
    else
        erparg1 = [];
        erparg2 = [];
    end;
    
    % adapt parameters
    % ----------------    
    if size(data,1) == 1 & size(data,3) == 1
        data = reshape(data, size(data,2)/length(sortvar), length(sortvar));
    end;
    if length(movewin) == 1, movewin(2) = 1; end;
    if length(decim)   == 1, decim(2)   = 0;   end;
    
    % sorting
    % -------
    if nargin > 2 & ~isempty(sortvar)
        [tmpsort tmpind] = sort(sortvar);
        sortind          = subjind(tmpind);
    else
        sortind          = subjind;
    end;
    
    % plot all erpimage and retreive values
    % -------------------------------------
    values = unique(sortind);
    nsubj = length(values);
    realdecim = decim(1);
    if ~plot_indiv, extra_args = { 'noshow', 'yes' };
    else            extra_args = { };
    end;
    for ival = 1:nsubj
        idxsubj = find(subjind == values(ival));
        p(ival) = length(idxsubj);
        
        if plot_indiv, figure; end;
        if movewin(1) >= 0
            if isempty(erparg1)
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                              times, titleim, movewin(1), -decim(1), varargin{:}, extra_args{:});
                %{ times, titleim, movewin(1), -decim(1), varargin{:} }
            else
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                              times, titleim, movewin(1), -decim(1), erparg1{3:end}, extra_args{:});
                %{  times, titleim, movewin(1), -decim(1), erparg1{3:end} }
            end;
        else
            % subtract moving average from erpimage()
            % ---------------------------------------
            if isempty(erparg1)
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, -movewin(1), -decim(1), varargin{:}, extra_args{:});
                [ tmpdata2 tmpvar2 ] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, 1, 0, varargin{:}, 'noshow', 'yes');
            else
                [ tmpdata tmpvar tmptrials tmplim axhndls,erp, ...
                  amps,cohers,cohsig,ampsig,outamps,phsangls,pshamp] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, -movewin(1), -decim(1), erparg1{3:end}, extra_args{:});
                [ tmpdata2 tmpvar2 ] = erpimage( data(:, idxsubj), sortvar( idxsubj), ...
                                                                  times, titleim, 1, 0, erparg1{3:end}, 'noshow', 'yes');
            end;
            if size(tmpdata2,2) > size(tmpdata,2)
                % randomly remove single trials
                % -----------------------------
                diflen   = size(tmpdata2,2) - size(tmpdata,2);
                indrm = find(shuffle([ ones(1, diflen) zeros(1, size(tmpdata,2)-diflen) ]));
                tmpdata2(:,indrm) = [];
                tmpvar2   (indrm) = [];
            end;
            maxlen = min(size(tmpdata,2), size(tmpdata2,2));
            tmpdata   = tmpdata2(:,1:maxlen)-tmpdata(:,1:maxlen);
            tmpvar    = tmpvar2(1:maxlen); % version non-smooth
            %tmpvar   = tmpvar(1:maxlen); % smooth version
            %phsangls  = phsangls(1:maxlen);
        end;

        % copy to variable regrouping all subject info
        % --------------------------------------------
        if ival == 1 
            realdecim  = size(tmpdata,2);
            outdata    = zeros(size(tmpdata,1), realdecim, nsubj);
            allsortvar = zeros(realdecim, nsubj);
            allphases  = zeros(realdecim, nsubj);
        end;
        
        if length(tmpdata) < size(outdata,2)
            outdata    = outdata   (:,1:length(tmpdata),:);
            allsortvar = allsortvar(1:length(tmpdata),:);
            if ~isempty(phsangls)
                allphases = allphases (1:length(tmpdata),:);
            end;
        end;
        minsize = min( size(outdata,2), size(tmpdata, 2) );
        outdata (:,1:minsize,ival) = tmpdata(:,1:minsize);     %/sum(sum(abs(tmpdata)));        
        allsortvar(1:minsize,ival) = tmpvar   (1:minsize)';
        
        if ~isempty(phsangls)
            tmp = movav(phsangls, [1:length(phsangls)], abs(movewin(1)), ...
                                    (length(phsangls)-abs(movewin(1))*2)/decim(1));
            allphases (:,ival) = tmp (1:size(allsortvar,1))';
        end;
    end;
    
    % uniform probability
    % -------------------
    uniform = p / sum(p);
    
    % resort trials and image them
    % ----------------------------
    newoutdata = mean(outdata,    3);
    newsortvar = mean(allsortvar, 2)';
    
    % plot global erpimage
    % --------------------
    timevect = linspace(tmplim(1), tmplim(2), size(outdata,1));
    if isempty(erparg2)
        % remove some arguments
        % ---------------------
        for index = length(varargin):-1:1
            if isstr(varargin{index}) & strcmpi(varargin{index}, 'align')
                varargin{index+1} = num2str(varargin{index+1});
            elseif isstr(varargin{index}) & strcmpi(varargin{index}, 'phasesort')
                varargin(index:index+1) = [];          
            elseif isstr(varargin{index}) & strcmpi(varargin{index}, 'plotamps')
                varargin(index) = [];          
            end;
        end;
        %{ times, titleim, movewin(2), -decim(2), varargin{:} }
        figure; erpimage( newoutdata, newsortvar, timevect, titleim, movewin(2), decim(2), 'nosort', varargin{:});
    else
        %{ times, titleim, movewin(2), -decim(2), erparg2{3:end} }
        figure; erpimage( newoutdata, newsortvar, timevect, titleim, movewin(2), decim(2), erparg2{3:end});
    end;
    return;
    
    % compute entropy
    % ---------------
    for sind = 1:nsubj % the loop below will select first
                       % the max trial among all subject
                       % second iteration will choose the second max...
                       % the loop has not been optimized
        
        % compute entropy array for this index
        % ------------------------------------
        subjcont = zeros(1, length(subjind));
        arrayp   = zeros(size(outdata,1), realdecim, nsubj);
        for index = 1:realdecim
            
            trials = outdata(:, index:realdecim:end); % nsubj trials (1 per subject)
            
            for indtime = 1:size(trials,1)
                [tmp indorder] = sort(trials(indtime,:));
                arrayp(indtime, index, indorder(sind)) = arrayp(indtime, index, indorder(1))+1;
            end;
            
        end;
        
        % compute entropy across time
        % ---------------------------
        arrayptime = squeeze(sum(arrayp,2));
        for indtime = 1:size(trials,1)
            p = arrayptime(indtime,:);
            p = p / sum(p);
            p(find(p == 0)) = [];
            subjamptime(indtime,sind) = -sum(p.*log(p)) / -sum(uniform.*log(uniform)) ;
        end;
        
        % compute entropy across trials
        % -----------------------------
        arrayptrial = squeeze(sum(arrayp,1));
        for indtrial = 1:size(arrayptrial,1)
            p = arrayptime(indtrial,:);
            p = p / sum(p);
            p(find(p == 0)) = [];
            subjamptrial(indtrial,sind) = -sum(p.*log(p)) / -sum(uniform.*log(uniform)) ;
        end;
    
        % global entropy
        % --------------
        p = squeeze(sum(sum(arrayp,1),2));
        p = p / sum(p);
        p(find(p == 0)) = [];
        globalent(sind) =  -sum(p.*log(p)) / -sum(uniform.*log(uniform));
    end;
    
    figure; plot(timevect, mean(subjamptime,2));
    figure; plot(1:realdecim, mean(subjamptrial,2));
    
    return
    
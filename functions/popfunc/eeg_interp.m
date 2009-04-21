% eeg_interp() - interpolate data channels
%
% Usage: EEGOUT = eeg_interp(EEG, badchans, method);
%
% Inputs: 
%     EEG      - EEGLAB dataset
%     badchans - [integer array] indices of channels to interpolate.
%                For instance, these channels might be bad.
%                [chanlocs structure] channel location structure containing
%                either locations of channels to interpolate or a full
%                channel structure (missing channels in the current 
%                dataset are interpolated).
%     method   - [string] griddata method used for interpolation
%                (default is 'invdist'). 'spherical' uses superfast
%                spherical interpolation. 'spacetime' uses griddata3 to
%                interpolate both in space and time (very slow).
% Output: 
%     EEGOUT   - data set with bad electrode data replaced by
%                interpolated data
%
% Author: Arnaud Delorme, CERCO, CNRS, Mai 2006-

% Copyright (C) Arnaud Delorme, CERCO, 2006, arno@salk.edu
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
% Revision 1.2  2008/04/16 17:34:45  arno
% added spherical and 3-D interpolation
%
% Revision 1.1  2006/09/12 18:46:30  arno
% Initial revision
%

function EEG = eeg_interp(ORIEEG, bad_elec, method)

    EEG = ORIEEG;
    if nargin < 2
        help eeg_interp;
        return;
    end;
    
    if nargin < 3
        method = 'spherical';
    end;

    if isstruct(bad_elec)
        
        % find missing channels
        % ---------------------
        if length(bad_elec) < length(EEG.chanlocs)
            bad_elec = [ EEG.chanlocs bad_elec ];
        end;
        if length(EEG.chanlocs) == length(bad_elec), return; end;
        
        lab1 = { bad_elec.labels };
        lab2 = { EEG.chanlocs.labels };
        [tmp badchans] = setdiff( lab1, lab2);
        fprintf('Interpolating %d channels...\n', length(badchans));
        goodchans      = sort(setdiff(1:length(bad_elec), badchans));
       
        % re-order good channels
        % ----------------------
        [tmp1 tmp2 neworder] = intersect( lab1, lab2 );
        [tmp1 ordertmp2] = sort(tmp2);
        neworder = neworder(ordertmp2);
        EEG.data = EEG.data(neworder, :, :);

        % looking at channels for ICA
        % ---------------------------
        %[tmp sorti] = sort(neworder);
        %{ EEG.chanlocs(EEG.icachansind).labels; bad_elec(goodchans(sorti(EEG.icachansind))).labels }
        
        % update EEG dataset (add blank channels)
        % ---------------------------------------
        if ~isempty(EEG.icasphere)
            
            [tmp sorti] = sort(neworder);
            EEG.icachansind = sorti(EEG.icachansind);
            EEG.icachansind = goodchans(EEG.icachansind);
            EEG.chaninfo.icachansind = EEG.icachansind;
            
            % TESTING SORTING
            %icachansind = [ 3 4 5 7 8]
            %data = round(rand(8,10)*10)
            %neworder = shuffle(1:8)
            %data2 = data(neworder,:)
            %icachansind2 = sorti(icachansind)
            %data(icachansind,:)
            %data2(icachansind2,:)
        end;
        % { EEG.chanlocs(neworder).labels; bad_elec(sort(goodchans)).labels }
        %tmpdata                  = zeros(length(bad_elec), size(EEG.data,2), size(EEG.data,3));
        %tmpdata(goodchans, :, :) = EEG.data;
        
        % looking at the data
        % -------------------
        %tmp1 = mattocell(EEG.data(sorti,1));
        %tmp2 = mattocell(tmpdata(goodchans,1));
        %{ EEG.chanlocs.labels; bad_elec(goodchans).labels; tmp1{:}; tmp2{:} }
        %EEG.data      = tmpdata;
        
        EEG.chanlocs  = bad_elec;

    else
        badchans  = bad_elec;
        goodchans = setdiff(1:EEG.nbchan, badchans);
    end;
    
    % find non-empty good channels
    % ----------------------------
    nonemptychans = find(~cellfun('isempty', { EEG.chanlocs.theta }));
    [tmp indgood ] = intersect(goodchans, nonemptychans);
    goodchans = goodchans( sort(indgood) );

    % scan data points
    % ----------------
    if strcmpi(method, 'spherical')
        % get theta, rad of electrodes
        % ----------------------------
        xelec = [ EEG.chanlocs(goodchans).X ];
        yelec = [ EEG.chanlocs(goodchans).Y ];
        zelec = [ EEG.chanlocs(goodchans).Z ];
        rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
        xelec = xelec./rad;
        yelec = yelec./rad;
        zelec = zelec./rad;
        xbad = [ EEG.chanlocs(badchans).X ];
        ybad = [ EEG.chanlocs(badchans).Y ];
        zbad = [ EEG.chanlocs(badchans).Z ];
        rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
        xbad = xbad./rad;
        ybad = ybad./rad;
        zbad = zbad./rad;
        
        EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
        %[tmp1 tmp2 tmp3 tmpchans] = spheric_spline_old( xelec, yelec, zelec, EEG.data(goodchans,1));
        %max(tmpchans(:,1)), std(tmpchans(:,1)), 
        %[tmp1 tmp2 tmp3 EEG.data(badchans,:)] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, EEG.data(goodchans,:));
        [tmp1 tmp2 tmp3 badchansdata] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, EEG.data);
        %max(EEG.data(goodchans,1)), std(EEG.data(goodchans,1))
        %max(EEG.data(badchans,1)), std(EEG.data(badchans,1))
        EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
    elseif strcmpi(method, 'spacetime') % 3D interpolation, works but x10 times slower
        pnts = size(EEG.data,2)*size(EEG.data,3);
        zgood = [1:pnts];
        zgood = repmat(zgood, [length(xgood) 1]);    zgood = reshape(zgood,prod(size(zgood)),1);
        xgood = repmat(xgood, [1 pnts]); xgood = reshape(xgood,prod(size(xgood)),1);
        ygood = repmat(ygood, [1 pnts]); ygood = reshape(ygood,prod(size(ygood)),1);
        tmpdata = reshape(EEG.data, prod(size(EEG.data)),1);
        zbad = 1:pnts;
        zbad = repmat(zbad, [length(xbad) 1]);     zbad = reshape(zbad,prod(size(zbad)),1);
        xbad = repmat(xbad, [1 pnts]); xbad = reshape(xbad,prod(size(xbad)),1);
        ybad = repmat(ybad, [1 pnts]); ybad = reshape(ybad,prod(size(ybad)),1);
        badchansdata = griddata3(ygood, xgood, zgood, tmpdata,...
                                              ybad, xbad, zbad, 'nearest'); % interpolate data                                            
    else 
        % get theta, rad of electrodes
        % ----------------------------
        [xbad ,ybad]  = pol2cart([EEG.chanlocs( badchans).theta],[EEG.chanlocs( badchans).radius]);
        [xgood,ygood] = pol2cart([EEG.chanlocs(goodchans).theta],[EEG.chanlocs(goodchans).radius]);

        fprintf('Points (/%d):', size(EEG.data,2)*size(EEG.data,3));
        badchansdata = zeros(length(badchans), size(EEG.data,2)*size(EEG.data,3));
        for t=1:(size(EEG.data,2)*size(EEG.data,3)) % scan data points
            if mod(t,100) == 0, fprintf('%d ', t); end;
            if mod(t,1000) == 0, fprintf('\n'); end;
            %for c = 1:length(badchans)
            %   [h EEG.data(badchans(c),t)]= topoplot(EEG.data(goodchans,t),EEG.chanlocs(goodchans),'noplot', ...
            %        [EEG.chanlocs( badchans(c)).radius EEG.chanlocs( badchans(c)).theta]);
            %end;
            tmpdata = reshape(EEG.data, size(EEG.data,1), size(EEG.data,2)*size(EEG.data,3) );
            [Xi,Yi,badchansdata(:,t)] = griddata(ygood, xgood , tmpdata(:,t)',...
                                                    ybad, xbad, method); % interpolate data                                            
        end
        fprintf('\n');
    end;
    
    tmpdata               = zeros(length(bad_elec), EEG.pnts, EEG.trials);
    tmpdata(goodchans, :,:) = EEG.data;
    tmpdata(badchans , :) = badchansdata;
    EEG.data = tmpdata;
    EEG.nbchan = size(EEG.data,1);
    EEG = eeg_checkset(EEG);
    
% -----------------
% spherical splines
% -----------------
function [x, y, z, Res] = spheric_spline_old( xelec, yelec, zelec, values);

SPHERERES = 20;
[x,y,z] = sphere(SPHERERES);
x(1:(length(x)-1)/2,:) = []; x = [ x(:)' ];
y(1:(length(y)-1)/2,:) = []; y = [ y(:)' ];
z(1:(length(z)-1)/2,:) = []; z = [ z(:)' ];

Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
Gsph  = computeg(x,y,z,xelec,yelec,zelec);

% equations are 
% Gelec*C + C0  = Potential (C unknow)
% Sum(c_i) = 0
% so 
%             [c_1]
%      *      [c_2]
%             [c_ ]
%    xelec    [c_n]
% [x x x x x]         [potential_1]
% [x x x x x]         [potential_ ]
% [x x x x x]       = [potential_ ]
% [x x x x x]         [potential_4]
% [1 1 1 1 1]         [0]

% compute solution for parameters C
% ---------------------------------
meanvalues = mean(values); 
values = values - meanvalues; % make mean zero
C = pinv([Gelec;ones(1,length(Gelec))]) * [values(:);0];

% apply results
% -------------
Res = zeros(1,size(Gsph,1));
for j = 1:size(Gsph,1)
    Res(j) = sum(C .* Gsph(j,:)');
end
Res = Res + meanvalues;
Res = reshape(Res, length(x(:)),1);

function [xbad, ybad, zbad, allres] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, values);

newchans = length(xbad);
numpoints = size(values,2);

%SPHERERES = 20;
%[x,y,z] = sphere(SPHERERES);
%x(1:(length(x)-1)/2,:) = []; xbad = [ x(:)'];
%y(1:(length(x)-1)/2,:) = []; ybad = [ y(:)'];
%z(1:(length(x)-1)/2,:) = []; zbad = [ z(:)'];

Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
Gsph  = computeg(xbad,ybad,zbad,xelec,yelec,zelec);

% compute solution for parameters C
% ---------------------------------
meanvalues = mean(values); 
values = values - repmat(meanvalues, [size(values,1) 1]); % make mean zero

values = [values;zeros(1,numpoints)];
C = pinv([Gelec;ones(1,length(Gelec))]) * values;
clear values;
allres = zeros(newchans, numpoints);

% apply results
% -------------
for j = 1:size(Gsph,1)
    allres(j,:) = sum(C .* repmat(Gsph(j,:)', [1 size(C,2)]));        
end
allres = allres + repmat(meanvalues, [size(allres,1) 1]);

% compute G function
% ------------------
function g = computeg(x,y,z,xelec,yelec,zelec)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +... 
                (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
                (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);

g = zeros(length(x(:)),length(xelec));
%dsafds
m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    L = legendre(n,EI);
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);    


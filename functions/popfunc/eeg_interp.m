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
%     method   - [string] method used for interpolation (default is 'spherical').
%                'invdist'/'v4' uses inverse distance on the scalp
%                'spherical' uses superfast spherical interpolation. 
%                'spacetime' uses griddata3 to interpolate both in space 
%                and time (very slow and cannot be interupted).
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

function EEG = eeg_interp(ORIEEG, bad_elec, method)

    if nargin < 2
        help eeg_interp;
        return;
    end;
    EEG = ORIEEG;
    
    if nargin < 3
        disp('Using spherical interpolation');
        method = 'spherical';
    end;

    % check channel structure
    tmplocs = ORIEEG.chanlocs;
    if isempty(tmplocs) || isempty([tmplocs.X])
        error('Interpolation require channel location');
    end;
    
    if isstruct(bad_elec)
        
        % add missing channels in interpolation structure
        % -----------------------------------------------
        lab1 = { bad_elec.labels };
        tmpchanlocs = EEG.chanlocs;
        lab2 = { tmpchanlocs.labels };
        [tmp tmpchan] = setdiff_bc( lab2, lab1);
        tmpchan = sort(tmpchan);
        if ~isempty(tmpchan)
            newchanlocs = [];
            fields = fieldnames(bad_elec);
            for index = 1:length(fields)
                if isfield(bad_elec, fields{index})
                    for cind = 1:length(tmpchan)
                        fieldval = getfield( EEG.chanlocs, { tmpchan(cind) },  fields{index});
                        newchanlocs = setfield(newchanlocs, { cind }, fields{index}, fieldval);
                    end;
                end;
            end;
            newchanlocs(end+1:end+length(bad_elec)) = bad_elec;
            bad_elec = newchanlocs;
        end;
        if length(EEG.chanlocs) == length(bad_elec), return; end;
        
        lab1 = { bad_elec.labels };
        tmpchanlocs = EEG.chanlocs;
        lab2 = { tmpchanlocs.labels };
        [tmp badchans] = setdiff_bc( lab1, lab2);
        fprintf('Interpolating %d channels...\n', length(badchans));
        if length(badchans) == 0, return; end;
        goodchans      = sort(setdiff(1:length(bad_elec), badchans));
       
        % re-order good channels
        % ----------------------
        [tmp1 tmp2 neworder] = intersect_bc( lab1, lab2 );
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
        goodchans = setdiff_bc(1:EEG.nbchan, badchans);
        oldelocs  = EEG.chanlocs;
        EEG       = pop_select(EEG, 'nochannel', badchans);
        EEG.chanlocs = oldelocs;
        disp('Interpolating missing channels...');
    end;

    % find non-empty good channels
    % ----------------------------
    origoodchans = goodchans;
    chanlocs     = EEG.chanlocs;
    nonemptychans = find(~cellfun('isempty', { chanlocs.theta }));
    [tmp indgood ] = intersect_bc(goodchans, nonemptychans);
    goodchans = goodchans( sort(indgood) );
    datachans = getdatachans(goodchans,badchans);
    badchans  = intersect_bc(badchans, nonemptychans);
    if isempty(badchans), return; end;
    
    % scan data points
    % ----------------
    if strcmpi(method, 'spherical')
        % get theta, rad of electrodes
        % ----------------------------
        tmpgoodlocs = EEG.chanlocs(goodchans);
        xelec = [ tmpgoodlocs.X ];
        yelec = [ tmpgoodlocs.Y ];
        zelec = [ tmpgoodlocs.Z ];
        rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
        xelec = xelec./rad;
        yelec = yelec./rad;
        zelec = zelec./rad;
        tmpbadlocs = EEG.chanlocs(badchans);
        xbad = [ tmpbadlocs.X ];
        ybad = [ tmpbadlocs.Y ];
        zbad = [ tmpbadlocs.Z ];
        rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
        xbad = xbad./rad;
        ybad = ybad./rad;
        zbad = zbad./rad;
        
        EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
        %[tmp1 tmp2 tmp3 tmpchans] = spheric_spline_old( xelec, yelec, zelec, EEG.data(goodchans,1));
        %max(tmpchans(:,1)), std(tmpchans(:,1)), 
        %[tmp1 tmp2 tmp3 EEG.data(badchans,:)] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, EEG.data(goodchans,:));
        [tmp1 tmp2 tmp3 badchansdata] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, EEG.data(datachans,:));
        %max(EEG.data(goodchans,1)), std(EEG.data(goodchans,1))
        %max(EEG.data(badchans,1)), std(EEG.data(badchans,1))
        EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
    elseif strcmpi(method, 'spacetime') % 3D interpolation, works but x10 times slower
        disp('Warning: if processing epoch data, epoch boundary are ignored...');
        disp('3-D interpolation, this can take a long (long) time...');
        tmpgoodlocs = EEG.chanlocs(goodchans);
        tmpbadlocs = EEG.chanlocs(badchans);
        [xbad ,ybad]  = pol2cart([tmpbadlocs.theta],[tmpbadlocs.radius]);
        [xgood,ygood] = pol2cart([tmpgoodlocs.theta],[tmpgoodlocs.radius]);
        pnts = size(EEG.data,2)*size(EEG.data,3);
        zgood = [1:pnts];
        zgood = repmat(zgood, [length(xgood) 1]);    
        zgood = reshape(zgood,prod(size(zgood)),1);
        xgood = repmat(xgood, [1 pnts]); xgood = reshape(xgood,prod(size(xgood)),1);
        ygood = repmat(ygood, [1 pnts]); ygood = reshape(ygood,prod(size(ygood)),1);
        tmpdata = reshape(EEG.data, prod(size(EEG.data)),1);
        zbad = 1:pnts;
        zbad = repmat(zbad, [length(xbad) 1]);     
        zbad = reshape(zbad,prod(size(zbad)),1);
        xbad = repmat(xbad, [1 pnts]); xbad = reshape(xbad,prod(size(xbad)),1);
        ybad = repmat(ybad, [1 pnts]); ybad = reshape(ybad,prod(size(ybad)),1);
        badchansdata = griddata3(ygood, xgood, zgood, tmpdata,...
                                              ybad, xbad, zbad, 'nearest'); % interpolate data                                            
    else 
        % get theta, rad of electrodes
        % ----------------------------
        tmpchanlocs = EEG.chanlocs;
        [xbad ,ybad]  = pol2cart([tmpchanlocs( badchans).theta],[tmpchanlocs( badchans).radius]);
        [xgood,ygood] = pol2cart([tmpchanlocs(goodchans).theta],[tmpchanlocs(goodchans).radius]);

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
            if strcmpi(method, 'invdist'), method = 'v4'; end;
            [Xi,Yi,badchansdata(:,t)] = griddata(ygood, xgood , double(tmpdata(datachans,t)'),...
                                                    ybad, xbad, method); % interpolate data                                            
        end
        fprintf('\n');
    end;
    
    tmpdata               = zeros(length(bad_elec), EEG.pnts, EEG.trials);
    tmpdata(origoodchans, :,:) = EEG.data;
    %if input data are epoched reshape badchansdata for Octave compatibility...
    if length(size(tmpdata))==3
        badchansdata = reshape(badchansdata,length(badchans),size(tmpdata,2),size(tmpdata,3));
    end
    tmpdata(badchans,:,:) = badchansdata;
    EEG.data = tmpdata;
    EEG.nbchan = size(EEG.data,1);
    EEG = eeg_checkset(EEG);

% get data channels
% -----------------
function datachans = getdatachans(goodchans, badchans);
      datachans = goodchans;
      badchans  = sort(badchans);
      for index = length(badchans):-1:1
          datachans(find(datachans > badchans(index))) = datachans(find(datachans > badchans(index)))-1;
      end;
        
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
    if ismatlab
        L = legendre(n,EI);
    else % Octave legendre function cannot process 2-D matrices
        for icol = 1:size(EI,2)
            tmpL = legendre(n,EI(:,icol));
            if icol == 1, L = zeros([ size(tmpL) size(EI,2)]); end;
            L(:,:,icol) = tmpL;
        end;
    end;
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);    


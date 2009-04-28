function [p] = eeg_colormap(p)

% eeg_colormap - generates EEG/ERP colormaps for eeg_toolbox
%
%   [p] = eeg_colormap(p)
%
%   Returns a colormap matrix [LENGTH, 3] into p.colorMap.map,
%   which can be used with matlab colormap command.  For example:
%
%  [p] = eeg_colormap; colormap(p.colorMap.map);
%
%   The p.colorMap field is a struct with the following fields:
%   
%   style:  'Gray'           grayscale (always linear)
%           'Red/Blue'       red/blue {default}
%           'Red/Blue/White' red/blue with white middle of width = WWIDTH
%   exp:    Exponent {1=linear, 2=squared [default],3=cubed,etc.}.
%           Higher exponents increase contrast between the extreme
%           values and the middle values.
%   length: size of colormap matrix {25 default}.
%           If not even, it will be rounded down to nearest even number
%   Cmin:   min RGB range, 0 <= RGB <=1; Cmin < Cmax {0 default}
%   Cmax:   max RGB range, 0 <= RGB <=1; Cmax > Cmin {1 default}
%   plot:   plot colormap, 1=plot, 0=no plot {0 default}
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2001, Darren.Weber_at_radiology.ucsf.edu, created
%           05/2001, Darren.Weber_at_radiology.ucsf.edu
%                    - added EXP argument and changed default from
%                      linear (EXP=1) to square (EXP=2).
%                    - ensured that Length is an even number
%                    - converted input/output to p struct
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),
    error('EEG_COLORMAP: No input p struct, use eeg_toolbox_defaults.');
else
    if ~isfield(p,'colorMap'), p.colorMap = []; end
    if ~isfield(p.colorMap,'style'),  p.colorMap.style   = 'Red/Blue/White'; end
    if ~isfield(p.colorMap,'exp'),    p.colorMap.exp     = 1;   end
	if ~isfield(p.colorMap,'Cmin'),   p.colorMap.Cmin    = 0;   end
	if ~isfield(p.colorMap,'Cmax'),   p.colorMap.Cmax    = 1;   end
	if ~isfield(p.colorMap,'plot'),   p.colorMap.plot    = 0;   end
end

% Ensure that Length is defined and an even number
if ~isfield(p.colorMap,'length'),
    p.colorMap.length  = 40;
else
    p.colorMap.length = fix(p.colorMap.length/2) * 2;
end


%A color map matrix may have any number of rows, but it must have
%exactly 3 columns.  Each row is interpreted as a color, with the
%first element specifying the intensity of red light, the second 
%green, and the third blue.  Color intensity can be specified on the
%interval 0.0 to 1.0.

switch p.colorMap.style
    
    case 'Gray'
        
        R = linspace(p.colorMap.Cmin,p.colorMap.Cmax,p.colorMap.length);
        Gmap = [R' R' R']; % R=G=B
        p.colorMap.map = Gmap;
        
    case 'Red/Blue'
        
        X = linspace(p.colorMap.Cmin,p.colorMap.Cmax,p.colorMap.length);
        R = X .^ p.colorMap.exp;
        G = R .* 0;
        B = fliplr(R);
        
        RBmap = [R' G' B'];
        p.colorMap.map = RBmap;
        
    case 'Red/Blue/White'
        
        Half = (p.colorMap.length / 2);
        
        % Calculate the required exponent min/max so that
        % min/max R1 below will be 0/Cmax-Cmin
        Expmin = 0;
        Expmax = (p.colorMap.Cmax - p.colorMap.Cmin) .^ ( 1 / p.colorMap.exp );
        
        % Generate an array of values from Expmin to Expmax, using
        % half the length of the colormap
        X  = linspace(Expmin,Expmax,Half);
        
        % Now raise it to the power of exp
        R1 = X .^ p.colorMap.exp;
        
        % Now add Cmin
        R1 = R1 + p.colorMap.Cmin;
        
        % Translate the curve origin to zero and
        % reflect it about the X axis
        tmp = 0 - (R1 - p.colorMap.Cmin);
        % Reflect it about the y axis
        tmp = fliplr(tmp);
        % Now add the Cmax value
        tmp = tmp + p.colorMap.Cmax;
        
        R1 = tmp(1:end-1);
        
        Cmax = ones(size(R1)) * p.colorMap.Cmax;
        
        R2 = fliplr(R1);
        
        R  = [ Cmax R2 ];   % Red array
        G  = [ R1 R2 ];     % Green array
        B  = fliplr(R);     % Blue array
        
        % Create colormap from Blue to Red
        Cmap = [B' G' R'];
        
        % Add the max value into the middle of the cmap to
        % even up the zero color point
        tmp1 = Cmap(1:Half-1,:);
        tmp2 = [p.colorMap.Cmax p.colorMap.Cmax p.colorMap.Cmax];
        tmp3 = Cmap(Half:end,:);
        Cmap = [ tmp1; tmp2; tmp3 ];
        
        p.colorMap.map = Cmap;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add pink or white tips to this colormap
        Expmin = 0;
        Expmax = (p.colorMap.Cmax - p.colorMap.Cmin) .^ ( 1 / 20 );
        
        X = linspace(Expmin,Expmax,length(Cmap));
        B = X .^ 10;
        R = fliplr(B);
        % for white tips
        %G = B + R; Addmap = [R' G' B'];
        % for pink/white tips
        %G = X .^ 40; G = G + fliplr(G); Addmap = [R' G' B'];
        % for pink tips
        %G = X .* 0; Addmap = [R' G' B'];
        % yellow, cyan tips
        G = B + R; Addmap = [B' G' R'];
        
        Addmap = (Cmap + Addmap);
        
        index = find(Addmap > p.colorMap.Cmax);
        Addmap(index) = p.colorMap.Cmax;
        
        p.colorMap.map = Addmap;
        
    case 'default',
        
        p.colorMap.map = jet;
        
    otherwise
        
        % get matlab colormap
        % this may not work under all circumstances, 
        % when Style is not a matlab colormap
        p.colorMap.map = colormap(p.colorMap.style);
        %R = linspace(Cmin,Cmax,Length);   B = fliplr(R);   G = zeros(size(R));
        %RBmap = [R' G' B'];
        %Map = RBmap;
end


% Plot colormap
if (p.colorMap.plot),
    figure('NumberTitle','off','Name','EEG Toolbox Colormap');
    colormap(p.colorMap.map); rgbplot(p.colorMap.map); colorbar
    axis tight
end

return

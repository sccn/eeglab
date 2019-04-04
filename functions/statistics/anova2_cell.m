%  anova2_cell() - compute F-values in cell array using ANOVA.
%
% Usage:
%    >> [FC FR FI dfc dfr dfi] = anova2_cell( data );
%
% Inputs:
%   data       = data consisting of PAIRED arrays to be compared. The last 
%                dimension of the data array is used to compute ANOVA.
% Outputs:
%   FC   - F-value for columns.
%   FR   - F-value for rows.
%   FI   - F-value for interaction.
%   dfc  - degree of freedom for columns.
%   dfr  - degree of freedom for rows.
%   dfi  - degree of freedom for interaction.
%
% Note: the advantage over the ANOVA2 function of Matlab statistical
%       toolbox is that this function works on arrays (see examples). Note
%       also that you still need the statistical toolbox to assess
%       significance using the fcdf() function. The other advantage is that
%       this function will work with complex numbers.
%
% Example:
%   a = { rand(1,10) rand(1,10) rand(1,10); rand(1,10) rand(1,10) rand(1,10) }
%   [FC FR FI dfc dfr dfi] = anova2_cell(a)
%   signifC = 1-fcdf(FC, dfc(1), dfc(2))
%   signifR = 1-fcdf(FR, dfr(1), dfr(2))
%   signifI = 1-fcdf(FI, dfi(1), dfi(2))
%
%   % for comparison 
%   anova2(  [ a{1,1}' a{1,2}' a{1,3}'; a{2,1}' a{2,2}' a{2,3}' ], 10) 
%
%   b = { [ a{1,1}; a{1,1} ] [ a{1,2}; a{1,2} ] [ a{1,3}; a{1,3} ];
%         [ a{2,1}; a{2,1} ] [ a{2,2}; a{2,2} ] [ a{2,3}; a{2,3} ] }
%   [FC FR FI dfc dfr dfi] = anova2_cell(b)
%
%   c{1,1} = reshape(repmat(b{1,1}, [2 1]),2,2,10);
%   c{1,2} = reshape(repmat(b{1,2}, [2 1]),2,2,10);
%   c{1,3} = reshape(repmat(b{1,3}, [2 1]),2,2,10);
%   c{2,3} = reshape(repmat(b{2,3}, [2 1]),2,2,10);
%   c{2,2} = reshape(repmat(b{2,2}, [2 1]),2,2,10);
%   c{2,1} = reshape(repmat(b{2,1}, [2 1]),2,2,10)
%   [FC FR FI dfc dfr dfi] = anova2_cell(c)
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005
%
% Reference:
%   Schaum's outlines in statistics (3rd edition). 1999. Mc Graw-Hill.

% Copyright (C) Arnaud Delorme
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [FC, FR, FI, freeC, freeR, freeI] = anova2_cell(data)
    
    % compute all means and all std
    % -----------------------------
    a = size(data,1);
    b = size(data,2);
    nd = myndims( data{1} );
    c  = size(data{1}, nd);
    sz = size(data{1});
    
    % dataabs if for complex data only
    % --------------------------------
    dataabs = data;
    if ~isreal(data{1})
        for i = 1:a
            for ii = 1:b
                dataabs{i,ii} = abs(data{i,ii});
            end
        end
    end
    
    VE = zeros( [ sz(1:end-1) 1], 'single' );
    m  = zeros( [ sz(1:end-1) size(data) ], 'single' );
    for i = 1:a
        for ii = 1:b
            tmpm = mymean(data{i,ii}, nd);
            switch nd
                case 1, m(i,ii) = tmpm;
                case 2, m(:,i,ii) = tmpm;
                case 3, m(:,:,i,ii) = tmpm;
                case 4, m(:,:,:,i,ii) = tmpm;
                case 5, m(:,:,:,:,i,ii) = tmpm;
                case 6, m(:,:,:,:,:,i,ii) = tmpm;
                case 7, m(:,:,:,:,:,:,i,ii) = tmpm;
                otherwise error('Dimension not supported');
            end
            VE = VE+sum( bsxfun(@minus, dataabs{i,ii}, tmpm).^2, nd);
        end
    end
    X  = mean(mean(m,nd+1),nd);
    Xj = mean(m,nd+1);
    Xk = mean(m,nd);
    VR = b*c*sum( bsxfun(@minus, Xj, X).^2, nd);
    VC = a*c*sum( bsxfun(@minus, Xk, X).^2, nd+1 );
    VI = c*sum( sum( bsxfun(@plus, bsxfun(@minus, bsxfun(@minus, m, Xj), Xk), X).^2, nd+1 ), nd );
        
% before bsxfun
%     VR = b*c*sum( (Xj-repmat(X, [ones(1,nd-1) size(Xj,nd  )])).^2, nd   );
%     VC = a*c*sum( (Xk-repmat(X, [ones(1,nd  ) size(Xk,nd+1)])).^2, nd+1 );
%     
%     Xj = repmat(Xj, [ones(1,nd)   size(m,nd+1) ]);
%     Xk = repmat(Xk, [ones(1,nd-1) size(m,nd)  1]);
%     VI = c*sum( sum( ( m - Xj - Xk + repmat(X, [ones(1,nd-1) size(m,nd) size(m,nd+1)]) ).^2, nd+1 ), nd );
    
%     if nd == 1
%         
%         VE = 0;
%         m  = zeros( size(data), 'single' );
%         for i = 1:a
%             for ii = 1:b
%                 m(i,ii) = mymean(data{i,ii});
%                 VE      = VE+sum( (dataabs{i,ii}-m(i,ii)).^2 );
%             end
%         end
%         X  = mean(mean(m));
%         Xj = mean(m,2);
%         Xk = mean(m,1);
%         VR = b*c*sum( (Xj-X).^2 );
%         VC = a*c*sum( (Xk-X).^2 );
%         
%         Xj = repmat(Xj, [1 size(m,2) ]);
%         Xk = repmat(Xk, [size(m,1)  1]);
%         VI = c*sum( sum( ( m - Xj - Xk + X ).^2 ) );
%         
%     elseif nd == 2
% 
%         VE = zeros( size(data{1},1),1, 'single');
%         m  = zeros( [ size(data{1},1) size(data) ], 'single' );
%         for i = 1:a
%             for ii = 1:b
%                 tmpm = mymean(data{i,ii}, 2);
%                 m(:,i,ii) = tmpm;
%                 VE        = VE+sum( (dataabs{i,ii}-repmat(tmpm, [1 size(data{i,ii},2)])).^2, 2);
%             end
%         end
%         X  = mean(mean(m,3),2);
%         Xj = mean(m,3);
%         Xk = mean(m,2);
%         VR = b*c*sum( (Xj-repmat(X, [1 size(Xj,2)])).^2, 2 );
%         VC = a*c*sum( (Xk-repmat(X, [1 1 size(Xk,3)])).^2, 3 );
%         
%         Xj = repmat(Xj, [1 1 size(m,3) ]);
%         Xk = repmat(Xk, [1 size(m,2)  1]);
%         VI = c*sum( sum( ( m - Xj - Xk + repmat(X, [1 size(m,2) size(m,3)]) ).^2, 3), 2 );
%         
%     elseif nd == 3
%         
%         VE = zeros( size(data{1},1), size(data{1},2), 'single' );
%         m  = zeros( [ size(data{1},1) size(data{1},2) size(data) ], 'single' );
%         for i = 1:a
%             for ii = 1:b
%                 tmpm = mymean(data{i,ii}, 3);
%                 m(:,:,i,ii) = tmpm;
%                 VE          = VE+sum( (dataabs{i,ii}-repmat(tmpm, [1 1 size(data{i,ii},3)])).^2, 3);
%             end
%         end
%         X  = mean(mean(m,4),3);
%         Xj = mean(m,4);
%         Xk = mean(m,3);
%         VR = b*c*sum( (Xj-repmat(X, [1 1 size(Xj,3)  ])).^2, 3 );
%         VC = a*c*sum( (Xk-repmat(X, [1 1 1 size(Xk,4)])).^2, 4 );
%         
%         Xj = repmat(Xj, [1 1 1 size(m,4) ]);
%         Xk = repmat(Xk, [1 1 size(m,3)  1]);
%         VI = c*sum( sum( ( m - Xj - Xk + repmat(X, [1 1 size(m,3) size(m,4)]) ).^2, 4 ), 3 );
%                 
%     else % nd == 4
%         
%         VE = zeros( size(data{1},1), size(data{1},2), size(data{1},3), 'single' );
%         m  = zeros( [ size(data{1},1) size(data{1},2) size(data{1},3) size(data) ], 'single' );
%         for i = 1:a
%             for ii = 1:b
%                 tmpm = mymean(data{i,ii}, 4);
%                 m(:,:,:,i,ii) = tmpm;
%                 VE            = VE+sum( (dataabs{i,ii}-repmat(tmpm, [1 1 1 size(data{i,ii},4)])).^2, 4);
%             end
%         end
%         X  = mean(mean(m,5),4);
%         Xj = mean(m,5);
%         Xk = mean(m,4);
%         VR = b*c*sum( (Xj-repmat(X, [1 1 1 size(Xj,4)  ])).^2, 4 );
%         VC = a*c*sum( (Xk-repmat(X, [1 1 1 1 size(Xk,5)])).^2, 5 );
%         
%         Xj = repmat(Xj, [1 1 1 1 size(m,5) ]);
%         Xk = repmat(Xk, [1 1 1 size(m,4)  1]);
%         VI = c*sum( sum( ( m - Xj - Xk + repmat(X, [1 1 1 size(m,4) size(m,5)]) ).^2, 5 ), 4 );
%                 
%     end
    
    SR2 = VR/(a-1);
    SC2 = VC/(b-1);
    SI2 = VI/(a-1)/(b-1);
    SE2 = VE/(a*b*(c-1));
    
    FR = SR2./SE2; % rows
    FC = SC2./SE2; % columns
    FI = SI2./SE2; % interaction
    
    freeR = [ a-1 a*b*(c-1) ];
    freeC = [ b-1 a*b*(c-1) ];
    freeI = [ (a-1)*(b-1) a*b*(c-1) ];

function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1,
            val = 2;
        elseif size(a,2) == 1,
            val = 1;
        else
            val = 2;
        end
    end; 
  
function res = mymean( data, varargin) % deal with complex numbers
    res = mean( data, varargin{:});
    if ~isreal(data)
        res = abs( res );
    end

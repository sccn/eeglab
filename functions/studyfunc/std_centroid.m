% std_centroid() - compute cluster centroid in EEGLAB dataset STUDY.
%                  Compute and store the centroid(s) (i.e., mean(s)) 
%                  for some combination of six measures on specified
%                  clusters in a STUDY. Possible measures include: scalp
%                  maps, ERPs, spectra, ERSPs, ITCs, dipole_locations
% Usage:    
%        >> [STUDY, centroid] = std_centroid(STUDY, ALLEEG, ...
%                                              clusters, measure1, measure2, ...);
%
% Inputs:
%   STUDY        - STUDY set 
%   ALLEEG       - ALLEEG dataset vector (else an EEG dataset) containing the STUDY
%                  datasets, typically created using load_ALLEEG().
%   clusters     - [vector] of cluster indices. Computes measure means for the 
%                  specified clusters. {deffault|[]: compute means for all 
%                  STUDY clusters} 
%   measure(s)   - ['erp'|'spec'|'scalp'|'dipole'|'itc'|'ersp'].   
%                  The measures(s) for which to calculate the cluster centroid(s):
%                     'erp'    ->  mean ERP of each cluster.
%                     'dipole' ->  mean dipole of each cluster.
%                     'spec'   ->  mean spectrum of each cluster (baseline removed).
%                     'scalp'  ->  mean topoplot scalp map of each cluster.
%                     'ersp'   ->  mean ERSP of each cluster. 
%                     'itc'    ->  mean ITC of each cluster. 
%                  If [], re-compute the centroid for whichever centroids 
%                  have previously been computed.
% Outputs:
%   STUDY        - input STUDY structure with computed centroids added. 
%                  If the requested centroids already exist, overwites them. 
%   centroid     - cell array of centroid structures, each cell corrasponding 
%                  to a different cluster requested in 'clusters' (above).
%                  fields of 'centroid' may include centroid.erp, centroid.dipole,
%                  etc. (as above). The structure is similar as the output
%                  of the std_readdata() function (with some fields
%                  about the cluster name and index missing).
% Examples:
%
%   >> [STUDY, centroid] = std_centroid(STUDY, ALLEEG,[], 'scalp'); 
%   % For each of the clusters in STUDY, compute a mean scalp map.
%   % The centroids are saved in the STUDY structure as entries in array
%   % STUDY.cluster(k).centroid.scalp. The centroids are also returned in 
%   % a cell array the size of the clusters (i.e., in: centroid(k).scalp).
%
%   >> [STUDY, centroid] = std_centroid(STUDY, ALLEEG,5,'spec','scalp'); 
%   % Same as above, but now compute only two centroids for Cluster 5. 
%   % The returned 'centroid' has two fields: centroid.scalp and centroid.spec
%
% Authors: Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, Feb 03, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Feb 03, 2005, hilit@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

% $Log: not supported by cvs2svn $
% Revision 1.26  2006/03/21 15:39:57  arno
% new .sets format
%
% Revision 1.25  2006/03/20 17:34:50  arno
% new std_clustread
%
% Revision 1.24  2006/03/14 02:55:44  scott
% help msg
%
% Revision 1.23  2006/03/12 04:29:24  arno
% G247
% transpose ERP
%
% Revision 1.22  2006/03/12 04:27:45  arno
% message
%
% Revision 1.21  2006/03/12 04:21:06  arno
% component polarity
%
% Revision 1.20  2006/03/11 07:04:18  arno
% header
%
% Revision 1.19  2006/03/11 00:29:15  arno
% header
%
% Revision 1.18  2006/03/10 18:22:17  arno
% renaming variables
%
% Revision 1.17  2006/03/10 16:27:08  arno
% new call for spectrum
%
% Revision 1.16  2006/03/10 16:00:00  arno
% reading ERP
% ./
%
% Revision 1.15  2006/03/10 15:56:04  arno
% changing log and var name
%
% Revision 1.14  2006/03/10 03:17:57  scott
% help msg -sm
%
% Revision 1.13  2006/03/09 23:48:43  arno
% fix new file format
%
% Revision 1.12  2006/03/09 23:31:29  arno
% log of frequencies
%
% Revision 1.11  2006/03/09 23:13:38  arno
% debug adding times field
%
% Revision 1.10  2006/03/09 23:11:00  arno
% adding field ersp_times to centroid
%
% Revision 1.9  2006/03/09 23:05:19  arno
% fixing call to std_readersp
%
% Revision 1.8  2006/03/09 22:26:25  arno
% reading ERSP to compute centroid
%
% Revision 1.7  2006/03/09 18:19:31  arno
% spectrum read
%
% Revision 1.6  2006/03/08 21:05:25  arno
% rename func
%
% Revision 1.5  2006/02/16 21:33:50  arno
% do not jump line for scalp centroid
%
% Revision 1.4  2006/02/15 22:12:42  arno
% revision
%

function [STUDY, centroid] = std_centroid(STUDY,ALLEEG, clsind, varargin);

 if nargin < 3
     help std_centroid;
     return
 end
 
 if isempty(clsind)
     for k = 2: length(STUDY.cluster) %don't include the ParentCluster
         if ~strncmpi('Notclust',STUDY.cluster(k).name,8) 
             % don't include 'Notclust' clusters
             clsind = [clsind k];
         end
     end
 end
 %defult values
erpC =0;
specC =0 ;
scalpC = 0;
dipoleC = 0;
itcC = 0;
erspC = 0;
 
commands = {};
if isempty(varargin)
    if isfield(STUDY.cluster(clsind(1)).centroid,'scalp')
        commands{end+1} = 'scalp';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'spec')
        commands{end+1} = 'spec';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'erp')
        commands{end+1} = 'erp';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'ersp')
        commands{end+1} = 'ersp';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'itc')
        commands{end+1} = 'itc';
    end
    if isfield(STUDY.cluster(clsind(1)).centroid,'dipole')
        commands{end+1} = 'dipole';
    end
else
    commands = varargin;
end

Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond = 1;
end 
centroid = cell(length(clsind),1);
fprintf('Computing ');
for k = 1:length(clsind)
    for l = 1:Ncond 
        for ind = 1:length(commands)
            ctr = commands{ind};
            switch ctr
                case 'scalp'
                    centroid{k}.scalp = 0; 
                    scalpC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('scalp ');
                    end
                case 'erp'
                    centroid{k}.erp{l} = 0; 
                    erpC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('erp ');
                    end
                case 'spec'
                    centroid{k}.spec{l} = 0; 
                    specC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('spectrum ');
                    end
                case 'ersp'
                    centroid{k}.ersp{l} = 0; 
                    centroid{k}.ersp_limits{l} = 0;
                    erspC =1;
                    if (l ==1) & (k ==1)
                        fprintf('ersp ');
                    end
                case 'itc'
                    centroid{k}.itc{l} = 0; 
                    centroid{k}.itc_limits{l} = 0;
                    itcC = 1;
                    if (l ==1) & (k ==1)
                        fprintf('itc ');
                    end                    
                case 'dipole'
                    dipoleC =1;
                    if (l ==1) & (k ==1)
                        fprintf('dipole ');
                    end
            end
        end
    end
end   
fprintf('centroid (only done once)\n');
if itcC | erspC | specC | erpC | scalpC
    for clust = 1:length(clsind) %go over all requested clusters
        for cond = 1:Ncond %compute for all conditions
            for k = 1:length(STUDY.cluster(clsind(clust)).comps) % go through all components
                comp  = STUDY.cluster(clsind(clust)).comps(k);
                abset = STUDY.cluster(clsind(clust)).sets(cond,k);
                if scalpC & cond == 1  %scalp centroid, does not depend on condition 
                    grid = std_readtopo(ALLEEG, abset, comp);
                    if isempty(grid)
                        return;
                    end
                    centroid{clust}.scalp = centroid{clust}.scalp + grid;
                end
                if erpC %erp centroid
                    [erp, t] = std_readerp(ALLEEG, abset, comp, STUDY.preclust.erpclusttimes);
                    fprintf('.');
                    if isempty(erp)
                        return;
                    end
                    if (cond==1) & (k==1)
                        all_erp = zeros(length(erp),length(STUDY.cluster(clsind(clust)).comps));
                    end
                    all_erp(:,k) = erp';
                    if k == length(STUDY.cluster(clsind(clust)).comps)
                        [all_erp pol] = std_comppol(all_erp);
                        centroid{clust}.erp{cond} = mean(all_erp,2);
                        centroid{clust}.erp_times = t;
                    end
                end
                if specC %spec centroid
                    [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
                    fprintf('.');
                    if isempty(spec)
                        return;
                    end
                    centroid{clust}.spec{cond} = centroid{clust}.spec{cond} + spec;
                    centroid{clust}.spec_freqs = f;
                end
                if erspC %ersp centroid
                    fprintf('.');
                    if cond == 1
                        tmpabset = STUDY.cluster(clsind(clust)).sets(:,k);
                        [ersp, logfreqs, timevals] = std_readersp(ALLEEG, tmpabset, comp, STUDY.preclust.erspclusttimes, ...
                                                                STUDY.preclust.erspclustfreqs );
                        if isempty(ersp)
                            return;
                        end
                        for m = 1:Ncond
                            centroid{clust}.ersp{m} = centroid{clust}.ersp{m} + ersp(:,:,m);
                            centroid{clust}.ersp_limits{m} = max(floor(max(max(abs(ersp(:,:,m))))), centroid{clust}.ersp_limits{m});
                        end
                        centroid{clust}.ersp_freqs  = logfreqs;
                        centroid{clust}.ersp_times = timevals;
                    end
                end
                if itcC %itc centroid
                    fprintf('.');
                    [itc, logfreqs, timevals] = std_readitc(ALLEEG, abset, comp, STUDY.preclust.erspclusttimes, ...
                                                                STUDY.preclust.erspclustfreqs );
                    if isempty(itc)
                        return;
                    end
                    centroid{clust}.itc{cond} = centroid{clust}.itc{cond} + itc;
                    centroid{clust}.itc_limits{cond} = max(floor(max(max(abs(itc)))), centroid{clust}.itc_limits{cond}); %ersp image limits 
                    centroid{clust}.itc_freqs  = logfreqs;
                    centroid{clust}.itc_times = timevals;
                end
            end
        end
        if ~scalpC
            fprintf('\n');
        end;
	end
end

if dipoleC %dipole centroid
    for clust = 1:length(clsind)
        max_r = 0;
        len = length(STUDY.cluster(clsind(clust)).comps);
        tmppos = 0;
        tmpmom = 0;
        tmprv = 0;
        ndip = 0;
        for k = 1:len 
            fprintf('.');
            comp  = STUDY.cluster(clsind(clust)).comps(k);
            abset = STUDY.cluster(clsind(clust)).sets(1,k);
            if ~isfield(ALLEEG(abset), 'dipfit')
               warndlg2(['No dipole information available in dataset ' num2str(abset) ], 'Aborting compute centroid dipole');
               return;
            end
            if ~isempty(ALLEEG(abset).dipfit.model(comp).posxyz)
                ndip = ndip +1;
                tmppos = tmppos + ALLEEG(abset).dipfit.model(comp).posxyz;
                tmpmom = tmpmom + ALLEEG(abset).dipfit.model(comp).momxyz;
                tmprv = tmprv + ALLEEG(abset).dipfit.model(comp).rv;
                if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                   if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                       load('-mat', ALLEEG(abset).dipfit.hdmfile);
                       max_r = max(max_r, max(vol.r));
                   else % old version of dipfit
                       max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
                   end
               end
            end
        end
        centroid{clust}.dipole.posxyz =  tmppos/ndip;
        centroid{clust}.dipole.momxyz =  tmpmom/ndip;
        centroid{clust}.dipole.rv =  tmprv/ndip;
        if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical') & (~isfield(ALLEEG(abset).dipfit, 'hdmfile')) %old dipfit
            centroid{clust}.dipole.maxr = max_r;
        end
        STUDY.cluster(clsind(clust)).centroid.dipole = centroid{clust}.dipole;
   end
end

%updat STUDY
for clust =  1:length(clsind) %go over all requested clusters
    for cond  = 1:Ncond
        ncomp = length(STUDY.cluster(clsind(clust)).comps);
        if scalpC & cond == 1%scalp centroid
            centroid{clust}.scalp  = centroid{clust}.scalp/ncomp;
            STUDY.cluster(clsind(clust)).centroid.scalp = centroid{clust}.scalp ;
        end
        if erpC
            STUDY.cluster(clsind(clust)).centroid.erp{cond} = centroid{clust}.erp{cond};
		    STUDY.cluster(clsind(clust)).centroid.erp_times = centroid{clust}.erp_times;
        end
		if specC
            centroid{clust}.spec{cond} = centroid{clust}.spec{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.spec{cond} = centroid{clust}.spec{cond};
		    STUDY.cluster(clsind(clust)).centroid.spec_freqs = centroid{clust}.spec_freqs;
        end
        if erspC %ersp centroid
            centroid{clust}.ersp{cond} = centroid{clust}.ersp{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.ersp{cond} = centroid{clust}.ersp{cond};
            STUDY.cluster(clsind(clust)).centroid.ersp_limits{cond} = floor(0.75*centroid{clust}.ersp_limits{cond}); 
            %[round(0.9*min(cell2mat({centroid{clust}.ersp_limits{cond,:}})))  round(0.9*max(cell2mat({centroid{clust}.ersp_limits{cond,:}})))];
            STUDY.cluster(clsind(clust)).centroid.ersp_freqs = centroid{clust}.ersp_freqs;
            STUDY.cluster(clsind(clust)).centroid.ersp_times = centroid{clust}.ersp_times;
        end
        if itcC
            centroid{clust}.itc{cond} = centroid{clust}.itc{cond}/ncomp;
            STUDY.cluster(clsind(clust)).centroid.itc{cond} = centroid{clust}.itc{cond} ;
            STUDY.cluster(clsind(clust)).centroid.itc_limits{cond} = floor(0.75*centroid{clust}.itc_limits{cond});%round(0.9*max(cell2mat({centroid{clust}.itc_limits{cond,:}})));
            STUDY.cluster(clsind(clust)).centroid.itc_freqs = centroid{clust}.itc_freqs;
            STUDY.cluster(clsind(clust)).centroid.itc_times = centroid{clust}.itc_times;
        end
        
    end
end
fprintf('\n');

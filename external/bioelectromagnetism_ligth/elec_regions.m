function [region] = elec_regions

% elec_regions - Defines a structure of electrode regions
% 
% Useage: [region] = elec_regions
% 
% The return struct is a 22 element array of
% structures with the fields:
% 
%   region(i).name  - string, region name/label
%   region(i).elec  - row vector of electrode numbers
%                     (eg, columns of ERP data)
%   region(i).label - cell array of electrode labels
% 
% Not all regions have the same number of electrodes.  The
% determination of regions was done at Flinders University
% with a Neuroscan/ECI electrode cap setup (circa 1996) to 
% correspond to scalp areas overlying various cortical 
% regions.  These regions were labeled as follows:
% 
%   IPF     - inferior prefrontal
%   IF      - inferior frontal
%   SPF     - superior prefrontal
%   SF      - superior frontal
%   SC      - superior central
%   IC      - inferior central
%   SP      - superior parietal
%   IP      - inferior parietal
%   AT      - anterior temporal
%   PT      - posterior temporal
%   OC      - occipital
% 
% The region names begin 'L_' or 'R_' indicating
% left or right hemisphere.  There are 22 regions in all.
% Note that midline sites are excluded from all regions.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% Created:  02/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

region( 1).name = 'L_IPF'; region( 1).elec = [1 4 6 7 121 122];          region( 1).labels = {'FP1','FP1p','AF7','AF3','aFP3','aFP1'};
region( 2).name = 'R_IPF'; region( 2).elec = [3 5 8 9 123 124];          region( 2).labels = {'FP2','FP2p','AF4','AF8','aFP2','aFP4'};

region( 3).name = 'L_IF';  region( 3).elec = [14 15 23 24 33];           region( 3).labels = {'F7','F5','F7p','F5p','FC5'};
region( 4).name = 'R_IF';  region( 4).elec = [21 22 29 30 39];           region( 4).labels = {'F6','F8','F6p','F8p','FC6'};

region( 5).name = 'L_SPF'; region( 5).elec = [10 11 16 17 26];           region( 5).labels = {'aF5','aF1','F3','F1','F1p'};
region( 6).name = 'R_SPF'; region( 6).elec = [12 13 19 20 27];           region( 6).labels = {'aF2','aF6','F2','F4','F2p'};

region( 7).name = 'L_SF';  region( 7).elec = [25 34 35 44 45];           region( 7).labels = {'F3p','FC3','FC1','aC3','aC1'};
region( 8).name = 'R_SF';  region( 8).elec = [28 37 38 46 47];           region( 8).labels = {'F4p','FC2','FC4','aC2','aC4'};

region( 9).name = 'L_SC';  region( 9).elec = [53 54 62 63 72];           region( 9).labels = {'C3','C1','C3p','C1p','CP1'};
region(10).name = 'R_SC';  region(10).elec = [55 56 64 65 74];           region(10).labels = {'C2','C4','C2p','C4p','CP2'};

region(11).name = 'L_IC';  region(11).elec = [42 43 52 60 61];           region(11).labels = {'aT7','aC5','C5','T7p','C5p'};
region(12).name = 'R_IC';  region(12).elec = [48 49 57 66 67];           region(12).labels = {'aC6','aT8','C6','C6p','T8p'};

region(13).name = 'L_SP';  region(13).elec = [71 81 82 91 99];           region(13).labels = {'CP3','aP3','aP1','P1','P1p'};
region(14).name = 'R_SP';  region(14).elec = [75 83 84 93 100];          region(14).labels = {'CP4','aP2','aP4','P2','P2p'};

region(15).name = 'L_IP';  region(15).elec = [70 79 80 89 90 98];        region(15).labels = {'CP5','aP7','aP5','P5','P3','P3p'};
region(16).name = 'R_IP';  region(16).elec = [76 85 86 94 95 101];       region(16).labels = {'CP6','aP6','aP8','P4','P6','P4p'};

region(17).name = 'L_AT';  region(17).elec = [31 32 50 51];              region(17).labels = {'FT9','FT7','T9','T7'};
region(18).name = 'R_AT';  region(18).elec = [40 41 58 59];              region(18).labels = {'FT8','FT10','T8','T10'};

region(19).name = 'L_PT';  region(19).elec = [68 69 87 88 102 111];      region(19).labels = {'TP9','TP7','P9','P7','TPO7','TO1'};
region(20).name = 'R_PT';  region(20).elec = [77 78 96 97 108 115];      region(20).labels = {'TP8','TP10','P8','P10','TPO8','TO2'};

region(21).name = 'L_OC';  region(21).elec = [103 104 109 112 116 118];  region(21).labels = {'PO7','PO3','aO1','O1','O1p','I1'};
region(22).name = 'R_OC';  region(22).elec = [106 107 110 114 117 120];  region(22).labels = {'PO4','PO8','aO2','O2','O2p','I2'};

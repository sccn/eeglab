function [scd] = eeg_lap2scd(lap,conductivity)

% eeg_lap2scd - scalp current density from Laplacian of scalp potential
%
% Useage: [scd] = eeg_lap2scd(lap,conductivity)
%
% where     'lap' is a scalp Laplacian matrix (ie, del^2(V)).
%           'conductivity' (> 0) is a scalp estimate. Commercial software 
%           (EMSE, SourceSignal; CURRY, Neuroscan) use 0.33 Siemens/m, 
%           which is the default here.  The conductivity units depend 
%           on the Laplacian units; V/m^2, uV/m^2, ?/m^2 is assumed.
%
% notes:    Current density is proportional to *minus* the Laplacian
%           of electric potential, so it has an opposite polarity to 
%           the Laplacian. This function employs the following equation:
%
%           div( current density ) = -1 ( conductivity * Laplacian(potential) )
%
%           div( current density ) = -1 * 0.33{S/m} * Laplacian {V/m^2}
%
%           The resulting values are Ampere/m^3 {A/m^3}.  If lap in uV/m^2, 
%           the result is in uA/m^3.  This result is often referred to as 
%           'scalp current density' or 'current source density' in the 
%           EEG literature.
%
% refs:     Geddes L, Baker L (1967) Med Biol Engr 5:271-293
%           Oostendorp T, van Oosterom A (1996) IEEE Trans Biomed Eng 43: 394-405
%           Perrin et al : http://www.lyon151.inserm.fr/unites/280_angl.html
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/2001, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist( 'conductivity', 'var' )
    conductivity = 0.33; % Siemens/meter
    fprintf ('Note: Scalp conductivity = %6.4f S/m; assume Laplacian in ?/m^2.\n', conductivity)
end

fprintf ('Scalp Current Density = -%6.4f * Laplacian(elec potential)\n', conductivity)
scd = -1 * conductivity * lap;

return

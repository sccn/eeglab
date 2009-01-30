%function [s,H]=fepi2gdf(fn);
% FEPI2GDF Freiburger epilepsy database converted into GDF format 

%	$Id: fepi2gdf.m,v 1.1 2009-01-30 06:04:41 arno Exp $
%	(C) 2005,2006,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.


if 1,
fid=fopen('epilepsydatabase.csv');
t=char(fread(fid,[1,inf],'uint8'));
fclose(fid);
[n,v,sa]=str2double(t,';',[10,13]);
end; 

fn = dir('*.asc');
x  = strvcat({fn.name});
id = unique(x(:,10:13),'rows');
for k=length(id)+(-2:0);%1:length(id)
        r=str2num(id(k,:));
        f = [fn(1).name(1:9),id(k,:)];

        sno = find(r==n(:,1)); 
        
        HDR.TYPE = 'GDF'; 
        HDR.FileName = [f,'.gdf'];
        HDR.Patient.Sex = strtok(sa{sno,2},'"'); 
        HDR.Patient.Age = n(sno,3); 
        HDR.T0 = clock; 
        HDR.Patient.Birthday = HDR.T0;
        HDR.Patient.Birthday(1) = HDR.T0(1) - HDR.Patient.Age;
        
        HDR.NS = 6; 
        HDR.DigMax  = (2^15-1)*ones(HDR.NS,1); 
        HDR.DigMin  = (-2^15)*ones(HDR.NS,1); 
        HDR.PhysMax = HDR.DigMax; 
        HDR.PhysMin = HDR.DigMin; 
        HDR.GDFTYP  = repmat(3,HDR.NS);
        HDR.FLAG.UCAL = 1; 
        HDR.SampleRate = 256; 
        HDR.Dur = 1/256;
        s = [];
        for k1 = 1:6, 
                t = load([f,'_',int2str(k1),'.asc']);
                s(:,k1) = t; 
        end;
        HDR.NRec = size(s,1); 
        HDR.SPR = 1; 
        
        HDR = sopen(HDR,'w'); 
        HDR = swrite(HDR,s); 
        HDR = sclose(HDR); 
        
        sview(HDR.FileName); 
        drawnow; 
end;

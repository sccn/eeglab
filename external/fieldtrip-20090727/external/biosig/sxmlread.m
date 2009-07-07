function [HDR]=sxmlopen(HDR,CHAN,arg4,arg5,arg6)
% SXMLOPEN reads XML files and tries to extract biosignal data
%
% This is an auxilary function to SOPEN. 
% Use SOPEN instead of SXMLREAD.
%

%
% HDR = sxmlopen(HDR);
%
% HDR contains the Headerinformation and internal data
%
% see also: SOPEN, SREAD, SSEEK, STELL, SCLOSE, SWRITE, SEOF


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%	$Id: sxmlread.m,v 1.1 2009-07-07 02:23:48 arno Exp $
%	(C) 2006 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


if strncmp(HDR.TYPE,'XML',3),
        if any(HDR.FILE.PERMISSION=='r'),
                fid = fopen(HDR.FileName,HDR.FILE.PERMISSION,'ieee-le');
                if strcmp(HDR.TYPE,'XML-UTF16'),
                        magic = char(fread(fid,1,'uint16'));
                        HDR.XML = char(fread(fid,[1,inf],'uint16'));
                elseif strcmp(HDR.TYPE,'XML-UTF8'),
                        HDR.XML = char(fread(fid,[1,inf],'char'));
                end;
                fclose(fid);
                HDR.FILE.FID = fid;
                if 1, try,
                        XML = xmltree(HDR.XML);
                        XML = convert(XML);
                        HDR.XML  =  XML; 
			HDR.TYPE = 'XML';
                catch
                        fprintf(HDR.FILE.stderr,'ERROR SOPEN (XML): XML-toolbox missing or invalid XML file.\n');
                        return;
                end;
                end;
		
		
                try,    % SierraECG  1.03  *.open.xml from PHILIPS
                        HDR.SampleRate = str2double(HDR.XML.dataacquisition.signalcharacteristics.samplingrate);
                        HDR.NS  = str2double(HDR.XML.dataacquisition.signalcharacteristics.numberchannelsvalid);
                        HDR.Cal = str2double(HDR.XML.reportinfo.reportgain.amplitudegain.overallgain);
                        HDR.PhysDim = 'uV';
                        HDR.Filter.HighPass = str2double(HDR.XML.reportinfo.reportbandwidth.highpassfiltersetting);
                        HDR.Filter.LowPass  = str2double(HDR.XML.reportinfo.reportbandwidth.lowpassfiltersetting);
                        HDR.Filter.Notch    = str2double(HDR.XML.reportinfo.reportbandwidth.notchfiltersetting);
                        
                        t = HDR.XML.reportinfo.reportformat.waveformformat.mainwaveformformat;
                        k = 0; 
                        HDR.Label=[];
                        while ~isempty(t),
                                [s,t] = strtok(t,' ');
                                k = k+1; 
                                HDR.Label{k, 1} = [s,' '];
                        end;
                        HDR.Patient.Id     = str2double(HDR.XML.patient.generalpatientdata.patientid);
                        tmp    = HDR.XML.patient.generalpatientdata.age;
                        if isfield(tmp,'years'),
	                        HDR.Patient.Age    = str2double(tmp.years);
			end
			if isfield(tmp,'dateofbirth')
				tmp = tmp.dateofbirth; 
				tmp(tmp=='-')=' ';
				HDR.Patient.Birthday([6,5,4]) = str2double(tmp);
			end;
                        
                        tmp    = HDR.XML.patient.generalpatientdata.sex;
                        HDR.Patient.Sex    = strncmpi(tmp,'Male',1) + strncmpi(tmp,'Female',1)*2; 
                        HDR.Patient.Weight = str2double(HDR.XML.patient.generalpatientdata.weight.kg);
                        HDR.Patient.Height = str2double(HDR.XML.patient.generalpatientdata.height.cm);
                        
                        HDR.VERSION = HDR.XML.documentinfo.documentversion;
                        HDR.TYPE = HDR.XML.documentinfo.documenttype;
                catch
                        
                try,    % FDA-XML Format
                        tmp   = HDR.XML.component.series.derivation;
                        if isfield(tmp,'Series');
                                tmp = tmp.Series.component.sequenceSet.component;
                        else    % Dovermed.CO.IL version of format
                                tmp = tmp.derivedSeries.component.sequenceSet.component;
                        end;
                        HDR.NS = length(tmp)-1;
                        HDR.NRec = 1; 
                        HDR.Cal = 1;
                        HDR.PhysDim = ' ';
                        HDR.SampleRate = 1;
                        HDR.TYPE = 'XML-FDA';     % that's an FDA XML file 
                catch
                
                try 	% GE Case8000 stress ECG

                	HDR.SampleRate = str2double(HDR.XML.StripData.SampleRate);
                	tmp = HDR.XML.ClinicalInfo.ObservationDateTime; 
                	HDR.T0 = [str2double(tmp.Year), str2double(tmp.Month), str2double(tmp.Day), str2double(tmp.Hour), str2double(tmp.Minute), str2double(tmp.Second)]; 

                	HDR.Patient.Id = HDR.XML.PatientInfo.PID; 
                	HDR.Patient.Name = 'X'; % [HDR.XML.PatientInfo.Name, ', ', HDR.XML.PatientInfo.GivenName];
                	HDR.Patient.Age = str2double(HDR.XML.PatientInfo.Age);
                	tmp = HDR.XML.PatientInfo.Gender; 
                	HDR.Patient.Sex = any(tmp(1)=='Mm') + any(tmp(1)=='Ff')*2;
                	HDR.Patient.Height = str2double(HDR.XML.PatientInfo.Height);
                	HDR.Patient.Weight = str2double(HDR.XML.PatientInfo.Weight); 
                	tmp = HDR.XML.PatientInfo.BirthDateTime; 
                	HDR.Patient.Birthday = [str2double(tmp.Year), str2double(tmp.Month), str2double(tmp.Day),0,0,0]; 
                	
                	tmp = HDR.XML.StripData.Strip; 
                	HDR.NS = length(tmp{1}.WaveformData);
               		tmax   = str2double(tmp{end}.Time.Minute)*60 + str2double(tmp{end}.Time.Second)+10;
                	data   = repmat(NaN,tmax*HDR.SampleRate,HDR.NS); 
                	for k  = 1:length(tmp);
	                	t = HDR.SampleRate*(str2double(tmp{k}.Time.Minute)*60 + str2double(tmp{k}.Time.Second));
	                	for k2 = 1:HDR.NS, 
		                	x = str2double(tmp{k}.WaveformData{k2}); 
		                	data(t+1:t+length(x),k2)=x(:); 
	                	end; 
                	end;
                	tmp = HDR.XML.ArrhythmiaData.Strip; 
                	for k  = 1:length(tmp);
	                	t = HDR.SampleRate*(str2double(tmp{k}.Time.Minute)*60 + str2double(tmp{k}.Time.Second));
	                	for k2 = 1:HDR.NS, 
		                	x = str2double(tmp{k}.WaveformData{k2}); 
		                	data(t+1:t+length(x),k2)=x(:); 
	                	end; 
                	end;
                	HDR.data = data - 2^12*(data>2^11); 
                	HDR.TYPE = 'native'; 
                	HDR.NRec = 1; 
                	HDR.SPR  = size(HDR.data,1); 
                	HDR.Calib = sparse(2:HDR.NS,1:HDR.NS,1); 
                	HDR.FLAG.UCAL = 1; 
                catch
                        fprintf(HDR.FILE.stderr,'Warning SOPEN (XML): File %s is not supported.\n',HDR.FileName);
                        return;
                end;
                end;
                end

                try
                	tmp=HDR.XML.componentOf.timepointEvent.componentOf.subjectAssignment.subject.trialSubject.subjectDemographicPerson.name; 
                	HDR.Patient.Name = sprintf('%s, %s',tmp.family, tmp.given);
		catch
                end
                
                
                
                HDR.Calib = sparse(2:HDR.NS+1,1:HDR.NS,HDR.Cal);
                HDR.FILE.OPEN = 1;
                HDR.FILE.POS  = 0;
        end;
end;

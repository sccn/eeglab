function [HDR,H1,h2]=opendicom(arg1,arg2,arg3,arg4,arg5,arg6)
% OPENDICOM is an auxillary function to SOPEN for 
% opening of DICOM files for reading ECG waveform data
% 
% Use SOPEN instead of OPENDICOM  
% 
% See also: fopen, SOPEN, 
%
% References: 
% [1] http://www.dclunie.com/dicom-status/status.html#BaseStandard2003
% [2] http://medical.nema.org/Dicom/supps/sup30_lb.pdf

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

%	$Revision: 1.1 $
%	$Id: opendicom.m,v 1.1 2009-01-30 06:04:41 arno Exp $
%	(C) 1997-2002, 2004 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

if nargin<2, 
        arg2='rb'; 
elseif ~any(arg2=='b');
        arg2= [arg2,'b']; % force binary open. 
end;
PERMISSION = arg2;

if isstruct(arg1) 
        HDR=arg1; 
        FILENAME=HDR.FileName;
else
        HDR.FileName=arg1;
        fprintf(2,'Warning OPENDICOM: the use of OPENDICOM is discouraged (OPENDICOM might disappear); please use SOPEN instead.\n');
end;


if any(PERMISSION=='r'),
	% Default Settings
	HDR.FLAG.implicite_VR = 1;
	HDR.Endianity = 'ieee-le';

	% Open file 
        HDR.FILE.FID = fopen(HDR.FileName,'r','ieee-le');
	id = fread(HDR.FILE.FID,132,'uchar');
	if ~all(id' == [zeros(1,128),abs('DICM')])
		status = fseek(HDR.FILE.FID,0,'bof');
	else
		HDR.FLAG.implicite_VR = 0;
        end;
		
        count = 0;
        [TAG,c] = fread(HDR.FILE.FID,2,'uint16');
        while ~feof(HDR.FILE.FID);
		TAG = [2^16,1]*TAG;
		if  HDR.FLAG.implicite_VR, 	% implicite VR
			LEN = fread(HDR.FILE.FID,1,'uint32');
		else  % Explicite VR
			[VR,c] = fread(HDR.FILE.FID,2,'uint16');
			ix = ['OB';'OW';'OF';'SQ';'UT';'UN']*[1;256];
			if ~c, 
				HDR.ERROR.status = -1;
				HDR.ERROR.message = sprintf('Error OPENDICOM: %s\n', HDR.FileName);
				return,
			end
			if ~any(VR(1)==ix)
				LEN = VR(2);
			else
				LEN = fread(HDR.FILE.FID,1,'uint32');
			end;	
		end;
		if (LEN==hex2dec('FFFFFFFF')),
			LEN = -1;
		end;

		%LEN = fread(HDR.FILE.FID,1,'uint32');
                count = count + 1;

		fprintf(1,'%03i\t%08x\t%04i\n',count,TAG,LEN);

%[count,LEN],dec2hex(TAG),		
		% read Value 
		if 0, 

		elseif (LEN>100000),  
			status = fseek(HDR.FILE.FID,LEN,0);
			if status, 
				fclose(HDR.FILE.FID);
				return;
			end;
			
		elseif (LEN<0), 
			fclose(HDR.FILE.FID);
			return;
		
		elseif (TAG==hex2dec('00000002')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.TAG00000002 = char(VAL');
		elseif (TAG==hex2dec('0000003a')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.TAG0000003a = VAL;

		elseif (TAG==hex2dec('00020000')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.GroupLength = VAL';
		elseif (TAG==hex2dec('00020001')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.FileMetaInformationVersion = VAL;
		elseif (TAG==hex2dec('00020002')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.UID = char(VAL');
		elseif (TAG==hex2dec('00020003')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.MediaStorage_SOP_InstanceUID = char(VAL)';
		elseif (TAG==hex2dec('00020010')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.TransferSyntaxUID = char(VAL');
		elseif (TAG==hex2dec('00020012')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.ImplementationClassUID = char(VAL');
		elseif (TAG==hex2dec('00020013')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uchar');
			HDR.DICOM.ImplementationVersionName = char(VAL');

		elseif (TAG==hex2dec('00080008')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ImageType = char(VAL);
		elseif (TAG==hex2dec('00080012')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.InstanceCreationDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080013')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.InstanceCreationTime = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080014')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.Instance_CreatorUID = char(VAL);
		elseif (TAG==hex2dec('00080016')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.SOP_ClassUID = char(VAL);
		elseif (TAG==hex2dec('00080018')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.SOP_InstanceUID = char(VAL);

		elseif (TAG==hex2dec('00080020')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.StudyDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080021')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.SeriesDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080022')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.AcquisitionDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080023')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ContentDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080024')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.OverlayDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080025')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.CurveDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('0008002A')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.AcquisitionDateTime = char(VAL(VAL~=abs('.')));

		elseif (TAG==hex2dec('00080030')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.StudyTime = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080031')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.SeriesTime = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080032')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.AcquisitionTime = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080033')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ContentTime = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080034')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.OverlayTime = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00080035')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.CurveTime = char(VAL(VAL~=abs('.')));

		elseif (TAG==hex2dec('00080050')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.AccessionNumber = char(VAL);
		elseif (TAG==hex2dec('00080060')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.Modality = char(VAL);
		elseif (TAG==hex2dec('00080070')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.Manufacturer = char(VAL);
		elseif (TAG==hex2dec('00080080')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.InstitutionName = char(VAL);
		elseif (TAG==hex2dec('00080090')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ReferringPhysiciansName = char(VAL);
		elseif (TAG==hex2dec('00081010')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.StationName = char(VAL);
		elseif (TAG==hex2dec('00081030')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.StudyDescription = char(VAL);
		elseif (TAG==hex2dec('00081040')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.InstitutionalDepartmentName = char(VAL);
		elseif (TAG==hex2dec('00081050')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.PerformingPhysiciansName = char(VAL);
		elseif (TAG==hex2dec('00081090')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ManufacturesModelName = char(VAL);
		elseif (TAG==hex2dec('0008114a')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ReferenceInstantSequence = char(VAL);
		elseif (TAG==hex2dec('0008114b')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ReferenceDescription = char(VAL);

		elseif (TAG==hex2dec('00100010')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.Patient.Name = char(VAL);
		elseif (TAG==hex2dec('00100020')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.PID = char(VAL);
		elseif (TAG==hex2dec('00100030')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.BirthDate = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00100030')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.BirthTime = char(VAL(VAL~=abs('.')));
		elseif (TAG==hex2dec('00100040')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.Patient.Sex = char(VAL);

		elseif (TAG==hex2dec('00181000')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.DeviceSerialNumber = char(VAL);
		elseif (TAG==hex2dec('00181020')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.SoftwareVersion = char(VAL);
		elseif (TAG==hex2dec('00181061')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.TriggerSource = char(VAL);
		elseif (TAG==hex2dec('00181067')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.ImageTriggerDelay = char(VAL);
		elseif (TAG==hex2dec('00181068')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.MultiplexGroupTimeOffset = char(VAL);
		elseif (TAG==hex2dec('00181069')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.TriggerTimeOffset = char(VAL);
		elseif (TAG==hex2dec('0018106a')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.SyncTrigger = char(VAL);
		elseif (TAG==hex2dec('0018106c')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.SyncChannel = char(VAL);
		elseif (TAG==hex2dec('0018106e')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.TriggerSamplePosition = char(VAL);
		elseif (TAG==hex2dec('00181800')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.AcqTimeSynchronized = char(VAL);
		elseif (TAG==hex2dec('00181801')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.TimeSource = char(VAL);
		elseif (TAG==hex2dec('00181802')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.TimeDistributionProtocol = char(VAL);
		elseif (TAG==hex2dec('00181810')),
			[VAL,c] = fread(HDR.FILE.FID,[1,LEN],'uchar');
			HDR.DICOM.AcqTimestamp = char(VAL);

		elseif (TAG==hex2dec('0020000d')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.StudyInstanceUID = char(VAL);
		elseif (TAG==hex2dec('0020000e')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.SeriesInstanceUID = char(VAL);
		elseif (TAG==hex2dec('00200010')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.StudyID = char(VAL);
		elseif (TAG==hex2dec('00200011')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.SeriesNumber = char(VAL);
		elseif (TAG==hex2dec('00200012')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.AcquisitionNumber = char(VAL);
		elseif (TAG==hex2dec('00200013')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.InstanceNumber = str2double(VAL);
		elseif (TAG==hex2dec('00200019')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.ItemNumber = char(VAL);
		elseif (TAG==hex2dec('00200054')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.TAG00200054 = VAL';

		elseif (TAG==hex2dec('00200200')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint8');
			HDR.DICOM.SychronizationFrame = char(VAL);

		elseif (TAG==hex2dec('00280010')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.Rows = VAL;
		elseif (TAG==hex2dec('00280011')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.Columns = VAL;
		elseif (TAG==hex2dec('00280012')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.Planes = VAL;
		elseif (TAG==hex2dec('00280030')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.PixelSpacing = VAL;
		elseif (TAG==hex2dec('00280100')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.BitsAllocated = VAL;
		elseif (TAG==hex2dec('00280101')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.BitsStored = VAL;
		elseif (TAG==hex2dec('00280102')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.HighBit = VAL;
		elseif (TAG==hex2dec('00280103')),
			[VAL,c] = fread(HDR.FILE.FID,1,'uint16');
			HDR.DICOM.PixelRepresentation = VAL;

		elseif (TAG==hex2dec('003a0004')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.WaveformOriginality = char(VAL); 
		elseif (TAG==hex2dec('003a0005')),	% waveform data
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.NS = VAL; 
		elseif (TAG==hex2dec('003a0010')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.SPR = VAL; 
		elseif (TAG==hex2dec('003a001a')),	% waveform data
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.SampleRate = VAL; 

		elseif (TAG==hex2dec('003a0020')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.MultiplexGroupLabel = VAL; 
		elseif (TAG==hex2dec('003a0200')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelDefinitionSequence = VAL; 
		elseif (TAG==hex2dec('003a0202')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelNumber = VAL; 
		elseif (TAG==hex2dec('003a0203')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelLabel = char(VAL); 
		elseif (TAG==hex2dec('003a0205')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelStatus = char(VAL); 
		elseif (TAG==hex2dec('003a0208')),	% waveform data
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelSource = VAL; 
		elseif (TAG==hex2dec('003a0209')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelSourceModifiersSequence = VAL; 
		elseif (TAG==hex2dec('003a020a')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.SourceWaveformSequence = VAL; 
		elseif (TAG==hex2dec('003a020c')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelDerivationDescription = VAL; 
		elseif (TAG==hex2dec('003a0210')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.Sensitivity = VAL; 
		elseif (TAG==hex2dec('003a0211')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.SensitivityUnits = VAL; 
		elseif (TAG==hex2dec('003a0212')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelSensitivityCorrectionFactor = char(VAL); 
		elseif (TAG==hex2dec('003a0213')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelBaseline = VAL; 
		elseif (TAG==hex2dec('003a0214')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelTimeSkew = VAL; 
		elseif (TAG==hex2dec('003a0215')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelSampleSkew = VAL; 
		elseif (TAG==hex2dec('003a0218')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ChannelOffset = VAL; 
		elseif (TAG==hex2dec('003a021a')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.SampleRate = VAL; 

		elseif (TAG==hex2dec('003a0220')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.Filter.LowFreq = VAL; 
		elseif (TAG==hex2dec('003a0221')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.Filter.HiFreq = VAL; 
		elseif (TAG==hex2dec('003a0222')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.Filter.NotchF = VAL; 
		elseif (TAG==hex2dec('003a0223')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.Filter.NotchB = VAL; 

		elseif (TAG==hex2dec('00400055')),
			[VAL,c] = fread(HDR.FILE.FID,LEN/2,'uint16');
			HDR.DICOM.AcquisitionContextModule = VAL;
		elseif (TAG==hex2dec('0040A043')),
			[VAL,c] = fread(HDR.FILE.FID,LEN/2,'uint16');
			HDR.DICOM.WaveformAnnotationModule = VAL;

		elseif (TAG==hex2dec('0040a0b0')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ReferencedWaveformChannel = VAL;
		elseif (TAG==hex2dec('0040a130')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.TemporalRangeType = VAL;
		elseif (TAG==hex2dec('0040a132')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ReferencedSamplePosition = VAL;
		elseif (TAG==hex2dec('0040a138')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ReferencedTimeOffset = VAL;
		elseif (TAG==hex2dec('0040a13a')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ReferencedDataTime = VAL;
		elseif (TAG==hex2dec('0040a180')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.AnnoationGroupName = char(VAL);
		elseif (TAG==hex2dec('0040a195')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.ConceptNameCodeSequenceModifier = VAL;
		elseif (TAG==hex2dec('0040B020')),
			[VAL,c] = fread(HDR.FILE.FID,LEN/2,'uint16');
			HDR.DICOM.WaveformAnnotationSequence = VAL;

		elseif (TAG==hex2dec('20330008'))
			%[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			[VAL,c] = fread(HDR.FILE.FID,LEN/2,'int16');
			HDR.DICOM.TAG20330008=VAL';
		elseif (TAG==hex2dec('42000000'))
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.TAG42000000=VAL';

		elseif (TAG==hex2dec('54000010')),	% waveform data
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			'54000010',
		elseif (TAG==hex2dec('54000100')),
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.WaveFormSequence = VAL; 
		elseif (TAG==hex2dec('54000110')),	% channel minimum
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.PhysMin = VAL;
		elseif (TAG==hex2dec('54000112')),	% channel maximum
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.PhysMax = VAL;
		elseif (TAG==hex2dec('54001004')),	% waveform bits allocation
			[HDR.bits,c] = fread(HDR.FILE.FID,LEN,'uint8');
		elseif (TAG==hex2dec('54001006')),	% waveform bits allocation
			[HDR.WaveformSampleInterpretation,c] = fread(HDR.FILE.FID,LEN,'uint8');
		elseif (TAG==hex2dec('54000100A')),	% padding
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.PaddingValue = VAL;
		elseif (TAG==hex2dec('54001010')),	% waveform data element
			[VAL,c] = fread(HDR.FILE.FID,LEN/2,'uint16');
			HDR.Data= VAL;

		elseif (TAG==hex2dec('7FE00010'))
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
			HDR.DICOM.TAG7fe00010=VAL';
		elseif (TAG==hex2dec('FFFEE000'))
			LEN = fread(HDR.FILE.FID,1,'uint32');
			if LEN~=hex2dec('FFFFFFFF')
				VAL = fread(HDR.FILE.FID,LEN,'uchar');
			else
				VAL = fread(HDR.FILE.FID,1,'uint32');
			end;
		else
			[VAL,c] = fread(HDR.FILE.FID,LEN,'uint8');
%			char(VAL'),
			fprintf(1,'ignored: TAG=%08x VAL=%s\n',TAG,char(VAL));
		end;	

                [TAG,c] = fread(HDR.FILE.FID,2,'uint16');
        end;
		
	
	if isfield(HDR.DICOM,'StudyDate') & isfield(HDR.DICOM,'StudyTime'),
		tmp = char(repmat(32,[1,20]));
		tmp([1:4,6,7,9,10]) = HDR.DICOM.StudyDate;
		tmp([12,13,15,16,18,19]) = HDR.DICOM.StudyTime(1:6);
		HDR.T0 = str2double(tmp);
	end;
	if isfield(HDR.DICOM,'Modality') & strncmp(HDR.DICOM.Modality,'ECG',3)
		tmp = char(repmat(32,[1,20]));
		tmp([1:4,6,7,9,10]) = HDR.DICOM.ContentDate;
		tmp([12,13,15,16,18,19]) = HDR.DICOM.ContentTime(1:6);
		HDR.T0 = str2double(tmp);
		HDR.TYPE = 'DICOM-ECG';
	end;

        fclose(HDR.FILE.FID);
end;        

% eeg_interp_sph_sline_test - script to test eeg_interp_sph_spline
%

clear

p = eeg_toolbox_defaults;

p.volt.path = 'e:\matlab\eeg_toolbox\eeg_example_data\';
p.volt.file = 'eeg_124ch_simulatedEMSE.txt'; % (400 rows x 125 columns)
p.volt.type = 'ascii';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load electrode co-ordinates

p = elec_open(p);
p.elec.plot = 1;

x = p.elec.data.x;
y = p.elec.data.y;
z = p.elec.data.z;

% % create spherical electrode positions
% r = 10
% fprintf('...generating spherical interpolation points\n');
% [x,y,z] = elec_sphere_points(16,24,r);
% 
% p.elec.data.x = x;
% p.elec.data.y = y;
% p.elec.data.z = z;
% 
% p.elec.data.Xsp = x;
% p.elec.data.Ysp = y;
% p.elec.data.Zsp = z;
% p.elec.data.Rsp = [r r r];
% 
% p.elec.n = 337;
%elec_plot(p)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create simulated potential data

p.volt.sampleHz = 1000;
p.volt.sampleMsec = 1;
p.volt.sampleTime = 1 * p.volt.sampleMsec;
p.volt.samplePoint = 1;
p.volt.epochStart = 0;
p.volt.epochEnd = 1;

p.volt.var = [];
p.volt.timeArray = 0:1:10;
p.volt.points = 10;
p.volt.channels = length(x);
p.volt.epochStart = 0;
p.volt.epochEnd = 10;
p.volt.sweeps = 1;
p.volt.peaks = [];

p.volt.data = repmat(p.elec.data.x',10,1);

V = p.volt.data(p.volt.samplePoint,:);

p.clickTimePoint = 0;
p = eeg_contours_engine(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spherical spline interpolation

FV = eeg_interp_sph_spline(V,[x y z]);

p.elec.n = length(FV.vertices);
p.elec.data.x = FV.vertices(:,1);
p.elec.data.y = FV.vertices(:,2);
p.elec.data.z = FV.vertices(:,3);
p.elec.data.Xsp = FV.vertices(:,1);
p.elec.data.Ysp = FV.vertices(:,2);
p.elec.data.Zsp = FV.vertices(:,3);

p.volt.data = repmat(FV.Cdata,10,1);

p.volt.file = 'spherical spline interpolation';

p = eeg_contours_engine(p);


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scd interpolation

FV = eeg_interp_sph_spline_scd(V,[x y z]);

p.elec.n = length(FV.vertices);
p.elec.data.x = FV.vertices(:,1);
p.elec.data.y = FV.vertices(:,2);
p.elec.data.z = FV.vertices(:,3);
p.elec.data.Xsp = FV.vertices(:,1);
p.elec.data.Ysp = FV.vertices(:,2);
p.elec.data.Zsp = FV.vertices(:,3);

p.volt.data = repmat(FV.Cdata,10,1);

p.volt.file = 'spherical spline scd interpolation';

p = eeg_contours_engine(p);

function elec_distance_bad(filename)

% elec_distance_bad - identify electrodes that were poorly digitised
%
% This utility is designed to identify electrodes that were
% poorly digitised by the Polhemus digitiser.  It calculates
% interelectrode distances and identifies electrodes with
% excessive distance (Mean +/- [2 * StDev]).
%
% The script creates output to text data files of nearest
% neighbour distances (*.nn) and outliers in both nearest
% neighbour distances and consecutive distances (*.bad).
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% History:  05/2000, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('filename','var')
        filename = 'eeg_example_data/elec_124_cart.txt';
    end

% Read electrode data, using utility 'elec_load'

    [elec,type,X,Y,Z] = elec_load(filename);
    
    % Find index values for electrodes, excluding PAN and centroid/ref
    electrodes = ones(size(type)) * 69;
    Index = find(type == electrodes);
    
    Ei = elec(Index);  Xi = X(Index);  Yi = Y(Index);  Zi = Z(Index);
    
% plot electrode position grid
    elec_meshplot(Xi,Yi,Zi);
    
    view(2);			% top view
%   view(-90, 0);		% left side view
%   view( 90, 0);		% right side view
%   view(-180,0);		% front view
%   view(   0,0);		% back view

% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate inter-electrode distances

  [dist_e, dist_d] = elec_distance(Ei,Xi,Yi,Zi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify distance between consecutive and nearest neighbour electrodes

  [C_e, C_d] = elec_distance_consec(dist_e, dist_d);
  [NN_e, NN_d] = elec_distance_nn(dist_e, dist_d);

  %clear dist_e, dist_d;

  NN_d = NN_d(:,2:end);  % remove column 1, all zero distances
  NN_e = NN_e(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for abnormal small or large distances

  fprintf('\n%s\n', 'Calculating outliers of electrode distances.');

  NN_bad_e = cell(0);
  NN_bad_d = [];

  NN_Mean = mean(NN_d(:,1));
  NN_StDev = std(NN_d(:,1));
  NN_low  = NN_Mean - (2 * NN_StDev);
  NN_high = NN_Mean + (2 * NN_StDev);

  for e = 1:length(NN_d(:,1))

	if(NN_d(e,1) > 0)

		if (NN_d(e,1) < NN_low)

			NN_bad_e(end + 1) = NN_e(e,1);
			NN_bad_d(end + 1) = NN_d(e,1);

		elseif (NN_d(e,1) > NN_high)

			NN_bad_e(end + 1) = NN_e(e,1);
			NN_bad_d(end + 1) = NN_d(e,1);
		end
	end
  end
  clear e;


  C_bad_e = cell(0);
  C_bad_d = [];

  C_Mean  = mean(C_d);
  C_StDev = std(C_d);
  C_low   = C_Mean - (2 * C_StDev);
  C_high  = C_Mean + (2 * C_StDev);

  for e = 1:length(C_d)

	if(C_d(e) > 0)

		if (C_d(e) < C_low)

			C_bad_e(end + 1) = C_e(e);
			C_bad_d(end + 1) = C_d(e);

		elseif (C_d(e) > C_high)

			C_bad_e(end + 1) = C_e(e);
			C_bad_d(end + 1) = C_d(e);
		end
	end
  end
  clear e;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% highlight bad electrodes on plot

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Nearest neighbour outliers


  if ~isempty(NN_bad_d)

	  Xb = [];  Yb = [];  Zb = [];

	  for n = 1:length(NN_bad_e)
	    e = str2num(char(NN_bad_e(n)));
	    index = strmatch(num2str(e(1)), char(elec), 'exact');
	    Xb(n) = X(index);
	    Yb(n) = Y(index);
	    Zb(n) = Z(index);
	  end
	  clear index e n;


	  plot3(Xb,Yb,Zb,'gx',Xb,Yb,Zb,'gd');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Consecutive electrode distance outliers

  if ~isempty(C_bad_d)

       Xb = [];  Yb = [];  Zb = [];

       for n = 1:length(C_bad_e)
         e = str2num(char(C_bad_e(n)));
         index = strmatch(num2str(e(1)), char(elec), 'exact');
         Xb(n) = X(index);
         Yb(n) = Y(index);
         Zb(n) = Z(index);
       end
       clear index e n;

       plot3(Xb,Yb,Zb,'bx',Xb,Yb,Zb,'bd');
  end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output 10 nearest neighbour electrode data

  [pathstr,name,ext,versn] = fileparts(filename);
  nnfile = fullfile(pathstr,strcat(name,'.nn'));

  FID = fopen(nnfile,'w');

  fprintf(  1,'%20s%10.4f\n','mean',NN_Mean,'stdev',NN_StDev,'low',NN_low,'high',NN_high);
  fprintf(FID,'%20s%10.4f\n','mean',NN_Mean,'stdev',NN_StDev,'low',NN_low,'high',NN_high);

  labels = str2num(char(NN_e(:,1)));
  nn = [labels NN_d(:,1)];

  for n = 2:10
	  labels = str2num(char(NN_e(:,n)));
	  nn(:,(end+1):(end+2)) = [labels(:,2) NN_d(:,n)];
  end

%  fprintf(  1,'%10d: %10d%10.4f%, 10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f\n',nn');
%  fprintf(FID,'%10d: %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f, %10d%10.4f\n',nn');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output bad electrode labels and distance values

  [pathstr,name,ext,versn] = fileparts(filename);
  bfile = fullfile(pathstr,strcat(name,'.bad'));

  FID = fopen(bfile,'w');

  if(NN_bad_d)

      NN_bad_en = str2num(char(NN_bad_e));
      NN_bad = [NN_bad_en(:,1) NN_bad_d']';
      fprintf('\n%s\n\n', 'Nearest neighbour summary stats:');
      fprintf('%-10s%10.4f\n','NN_Mean',NN_Mean,'NN_StDev',NN_StDev,'NN_low',NN_low,'NN_high',NN_high);
      fprintf('\n%s\n\n', 'Nearest neighbour outliers (green):');
      fprintf('%10d%10.4f\n', NN_bad);

      fprintf(FID,'\n%s\n\n', 'Nearest neighbour summary stats:');
      fprintf(FID,'%-10s%10.4f\n','NN_Mean',NN_Mean,'NN_StDev',NN_StDev,'NN_low',NN_low,'NN_high',NN_high);

      fprintf(FID,'\n%s\n\n', 'Nearest neighbour outliers:');
      fprintf(FID,'%10d%10.4f\n', NN_bad);


  end

  if(C_bad_d)

      C_bad_en = str2num(char(C_bad_e));
      C_bad = [C_bad_en(:,1) C_bad_d'];
      fprintf('\n%s\n\n', 'Consecutive electrode summary stats:');
      fprintf('%-10s%10.4f\n','C_Mean',C_Mean,'C_StDev',C_StDev,'C_low',C_low,'C_high',C_high);
      fprintf('\n%s\n\n', 'Consecutive electrode outliers (blue):');
      fprintf('%10.0f%10.4f\n', C_bad);

      fprintf(FID,'\n%s\n\n', 'Consecutive electrode summary stats:');
      fprintf(FID,'%-10s%10.4f\n','C_Mean',C_Mean,'C_StDev',C_StDev,'C_low',C_low,'C_high',C_high);

      fprintf(FID,'\n%s\n\n', 'Consecutive electrode outliers:');
      fprintf(FID,'%10d%10.4f\n', C_bad);

  end

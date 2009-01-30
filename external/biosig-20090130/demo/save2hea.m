[s0,HDR]=sload('cardiocontrol.com000001.scp');
H2 = HDR; 

H2.TYPE   = 'MIT';
H2.GDFTYP = repmat(3,1,HDR.NS); 
R = y2res(s0);

H2.Off     = zeros(1,HDR.NS); 
H2.Cal     = R.QUANT;
H2.PhysMin = R.Min;
H2.PhysMax = R.Max;
H2.DigMin  = H2.PhysMin./R.QUANT; 
H2.DigMax  = H2.PhysMax./R.QUANT; 

H2 = sopen(H2,'w'); 
H2 = swrite(H2,s0);
H2 = sclose(H2); 

[s3,H3] = sload(H2.FileName); 



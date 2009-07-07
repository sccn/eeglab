function crc16 = crc16eval(D)
% CRC16EVAL cyclic redundancy check with the polynomiaL x^16+x^12+x^5+1  
% i.e. CRC-CCITT http://en.wikipedia.org/wiki/Crc16 

	D = uint16(D);

	crchi = 255;
	crclo = 255;

	t = '00102030405060708191a1b1c1d1e1f112023222524272629383b3a3d3c3f3e32434041464744454a5b58595e5f5c5d53626160676665646b7a79787f7e7d7c74858687808182838c9d9e9f98999a9b95a4a7a6a1a0a3a2adbcbfbeb9b8bbbab6c7c4c5c2c3c0c1cedfdcdddadbd8d9d7e6e5e4e3e2e1e0effefdfcfbfaf9f8f9181b1a1d1c1f1e110003020504070608393a3b3c3d3e3f30212223242526272b5a59585f5e5d5c53424140474645444a7b78797e7f7c7d72636061666764656d9c9f9e99989b9a95848786818083828cbdbebfb8b9babbb4a5a6a7a0a1a2a3afdedddcdbdad9d8d7c6c5c4c3c2c1c0cefffcfdfafbf8f9f6e7e4e5e2e3e0e1e';
	crc16htab = hex2dec(reshape(t,2,length(t)/2)');

	t = '0021426384a5c6e708294a6b8cadceef31107352b594f7d639187b5abd9cffde62432001e6c7a4856a4b2809eecfac8d53721130d7f695b45b7a1938dffe9dbcc4e586a740610223cced8eaf48690a2bf5d4b79671503312fddcbf9e79583b1aa687e4c522036041ae8feccd2a0b684997b6d5f4133251709fbeddfc1b3a597888a9caeb0c2d4e6f80a1c2e304254667b998fbda3d1c7f5eb190f3d235147756eacba8896e4f2c0de2c3a08166472405dbfa99b85f7e1d3cd3f291b0577615344c6d0e2fc8e98aab44650627c0e182a37d5c3f1ef9d8bb9a75543716f1d0b3922e0f6c4daa8be8c926076445a283e0c11f3e5d7c9bbad9f81736557493b2d1f0';
	crc16ltab = hex2dec(reshape(t,2,length(t)/2)');

	for k = 1:length(D),
		ix = bitxor(crchi,D(k))+1;
		crchi = bitxor(crclo,crc16htab(ix));
		crclo = crc16ltab(ix);
	end;
	crc16 = double(crchi)*256+crclo;
	
end;	%%%%% crc16eval %%%%%



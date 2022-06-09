function genomat = PlinkRead_binary2(nsubj,snps,fileprefix,subjs)
persistent geno_values

if ~exist('subjs','var'), subjs = 1:nsubj; end % By defaults, include all subjects

% Written by Chun 2015
if ~issorted(snps), error('PlinkRead_binary2 expect sorted list of snps'); end;
nsnp = length(snps);

% bit shift to generate the genovalue matrix
bedprefix = sprintf('%s.bed', fileprefix);

if isempty(geno_values)
    geno_values = zeros(256,4,'int8');
    geno_code = [-1,1,0,2];
    shiftind = [0,2,4,6];
    indvec=zeros(1,4);

    for i = 1:256;
        ishift = int16(i-1);
        for j = 1:4;
            indvec(j) = bitand(bitsra(ishift,shiftind(j)),3) ;
        end
        indvec(indvec == 0) = 4;
        geno_values(i,:) = geno_code(indvec);
    end
end

% Read in the binary file
bedid = fopen(bedprefix);
genobin = uint16(fread(bedid, 3));

% Check magic number
if genobin(1) ~= 108;
	error('- Not a valid Plink BED file \r\n');
elseif genobin(2) ~= 27;
	error('- Not a valid Plink BED file \r\n');
elseif genobin(3) ~= 1;
	error('- Not in SNP-major format \r\n');
end

n_bytes_per_snp = ceil(nsubj/4);
byte0 = (subjs(1)-1)/4; % Make sure subjs(1) is integer multiple of 4
byte1 = ceil((subjs(end))/4);
n_bytes_per_read = byte1-byte0;
genomat = zeros(length(subjs),nsnp,'int8');
for i = 1:nsnp
    if mod(i,1e6)==0
      fprintf(1,'snpi=%d/%d (%s)\n',i,nsnp,datestr(now));
    end
    fseek(bedid, 3 + (snps(i) - 1) * n_bytes_per_snp+byte0, 'bof');
    genobin = uint16(fread(bedid, n_bytes_per_read));
    if length(genobin) ~= n_bytes_per_read, error('-- Invalid number of entries from the bed \r\n'); end
    tmp_values = geno_values(genobin + 1, :)';
    tmp_values = tmp_values(:);
    genomat(:,i) = tmp_values(1:length(subjs));
end
fclose(bedid);


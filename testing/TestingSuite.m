addpath('/data/Linderman/FastPCA2/Matlab/')
addpath('/data/Linderman/software/software/FastPcaMatlabTygert/')

errs = []
%%%%%%%%%%%%%%%%%%%%%%%%%EIGENSTRATGENO: with imputation
inputFile = 'example_with_impute.eigenstratgeno';
k = 5;
l = 5;
m = 7;
n = 5;

%%%%%%%%%%Use C 
cmd = sprintf('../bin/fastpca_ooc.xx -k %d -eigenInCoreInput %s -binaryOutput example.m7n5k5_incore.binmatrix -eigenOutput example.m7n5k5_incore.eigenout -binA centeredMatrix_incore.binmatrix -d -c -mem 1.4e+2 -l %d', k, inputFile,  l)
system(cmd);

disp(cmd);
cmd = sprintf('../bin/fastpca_ooc.xx -k %d -eigenInput %s -binaryOutput example.m7n5k5.binmatrix -eigenOutput example.m7n5k5.eigenout -binA centeredMatrix.binmatrix -d -c -mem 1.4e+2 -l %d', k, inputFile,  l)




system(cmd);

%%%%%%%%%%Use Matlab
fileID=fopen(inputFile);
frewind(fileID);
formatString = repmat('%c',1,n);
fileText = textscan(fileID,formatString,'CollectOutput',1); 
A= cast(fileText{1,1}(:,:),'double')-48;

%Imputation of columns, so it can work with the matlab pca
for i=1:m,
        I = find (A(i,:)==9 );
        NonMissing = setdiff(1:n,I);
        meanRow = mean(A(i, NonMissing));
        A(i,I) = meanRow;
end
Ac = A - mean(A,2) * ones(1,n) ;
[U_g,S_g,V_g] = loadFastPcaBinaries('example.m7n5k5_incore.binmatrix',7,5,5);
normImputedSmall = diffsnorm(Ac, U_g,S_g,V_g);
errs = [errs normImputedSmall] 
if (normImputedSmall > 10E-14)
	error('Test Failed: Small matrix with imputation incore');
end

[U_g,S_g,V_g] = loadFastPcaBinaries('example.m7n5k5.binmatrix',7,5,5);
normImputedSmall = diffsnorm(Ac, U_g,S_g,V_g);
errs = [errs normImputedSmall] 
if (normImputedSmall > 10E-14)
	error('Test Failed: Small matrix with imputation out of core');
end
%Can also compare U_g with U_m below, within a factor of -1 should be the same
%[U_m, S_m, V_m] = pca(Ac,k, true);
%diffsnorm(Ac,U_g,S_g,V_g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Binary input, using an artificial matrix of rank 20
m=1000;
n=3000;
k=20;
binOutput= 'example.m1000n3000k20.binmatrix';
binInput = '/data/BigDataSandbox/1000x3000.binmatrix'
for memory=[1e+6 1e+5 30e+4 1e+9]
%%Check the -binaryInput

	cmd = sprintf('../bin/fastpca_ooc.xx -k %d -binaryInput %s  -binaryOutput  %s -d -c -mem %d -m %d -n %d',k, binInput, binOutput, memory, m,n);
	disp (cmd)
	system(cmd);

	A = loadFastPcaBinary(binInput, m,n);
	Ac = A - mean(A,2) * ones(1,n) ;
	[U_g,S_g,V_g] = loadFastPcaBinaries(binOutput,m,n, k);
	normBinaryArtificial = diffsnorm(Ac,U_g,S_g,V_g);
	errs = [errs normBinaryArtificial] 
	if (normBinaryArtificial > 10E-10)
		error('Test Failed: Artificial Binary matrix');
	end

end


%%%%Test
k=20
cmd = sprintf('../bin/fastpca_ooc.xx -k %d -binaryInput %s  -binaryOutput  %s -d -c -mem %d -m %d -n %d',k, binInput, binOutput, memory, m,n);
disp (cmd)
system(cmd);

A = loadFastPcaBinary(binInput, m,n);
Ac = A - mean(A,2) * ones(1,n) ;
[U_t,S_t,V_t] = pca(Ac, k, true);
normBinaryArtificialT = diffsnorm(Ac,U_t,S_t,V_t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BED Input: Did not have a Matlab BED reader, so used convertf to convert it into an eigenstratgeno
%%%%Check the imputation and centering
fileID=fopen('/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_1000_62318/eur_1000_62318.2.eigenstratgeno');
bedFiles = '/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_1000_62318/eur_1000_62318'
binOutput= 'eur_1000_62318_k20.9.bim.binmatrix'
m = 958; %SNPs
n =  62318; %Patients
k = 30;

frewind(fileID);
formatString = repmat('%c',1,n);
EUR_txt = textscan(fileID,formatString,'CollectOutput',1);
A= cast(EUR_txt{1,1}(:,:),'double')-48;

for i=1:m,
        I = find (A(i,:)==9 );
        NonMissing = setdiff(1:n,I);
        meanRow = mean(A(i, NonMissing));
        A(i,I) = meanRow;
end
Ac = A - mean(A,2) * ones(1,n) ;



%%%%Check the PCA for a large BED dataset
bedFiles = '/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_20000_62318/eur_20000_62318'
binOutput= 'eur_20000_62318_k20.bim.binmatrix'
binA = 'centeredMatrix.20000_62318.binmatrix'
m = 19826; %SNPs
n =  62318; %Patients
k = 30;
memory = 30e+9

%Check the PCA

for memory=[15e+9 5e+9  30e+9 ]
	for k =[10 20 40]
		cmd = sprintf('../bin/fastpca_ooc.xx -k %d -bedI  %s -binaryOutput  %s -binA %s -mem %d -c ', k, bedFiles, binOutput, binA, memory);
		disp (cmd)
		system(cmd);

		Ag = loadFastPcaBinary(binA, m,n);
		[U_g,S_g,V_g] = loadFastPcaBinaries(binOutput,m, n, k);
		norm_g = diffsnorm(Ag,U_g,S_g,V_g)

		[U_t,S_t,V_t] = pca(Ag, k,true);
		norm_t = diffsnorm(Ag,U_t,S_t,V_t)


		if (abs(norm_t - norm_g) > 30) 
			error('Matlab and PCA norms are quite different');
		end
	end
end

if (norm(abs(Ag) - abs(Ac)> 10e-10) )
	error ('The imputation and centering is inaccurate')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSV Output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSV Output

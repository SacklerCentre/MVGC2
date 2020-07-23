
if ~exist('n',    'var'), n     = 11;   end
if ~exist('r',    'var'), r     = 15;   end
if ~exist('rhoa', 'var'), rhoa  = 0.9;  end
if ~exist('rmi',  'var'), rmi   = 0.5;  end
if ~exist('fres', 'var'), fres  = [];   end
if ~exist('tol',  'var'), tol   = [];   end
if ~exist('maxi', 'var'), maxi  = [];   end
if ~exist('seed', 'var'), seed  = 0;    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng_seed(seed);

[A,C,K] = iss_rand(n,r,rhoa);
V = corr_rand(n,rmi);
ssinfo = ss_info(A,C,K,V);

if isempty(fres)
	maxfres = 8096;
    fres = 2^nextpow2(ssinfo.acdec);
	if fres > maxfres % adjust to taste
		fprintf(2,'ARNING: large frequency resolution = %d. Resetting to %d\n\n',fres,maxfres);
		fres = maxfres;
	else
		fprintf('Using frequency resolution %d\n\n',fres);
	end
end

S = ss_to_cpsd(A,C,K,V,fres);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1 = ss_to_pwcgc(A,C,K,V);

F2 = cpsd_to_pwcgc(S,tol,maxi);

fprintf('F1 =\n'); disp(F1);
fprintf('F2 =\n'); disp(F2);
RE = 2*maxabs(F1(:)-F2(:))/(maxabs(F1(:))+maxabs(F2(:)));
fprintf('RE = %e\n\n',RE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = ss_to_spwcgc(A,C,K,V,fres);

f2 = cpsd_to_spwcgc(S,tol,maxi);

re = 2*maxabs(f1(:)-f2(:))/(maxabs(f1(:))+maxabs(f2(:)));

fprintf('re = %e\n\n',re);

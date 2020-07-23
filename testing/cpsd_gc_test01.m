
if ~exist('n',    'var'), n     = 11;   end
if ~exist('r',    'var'), r     = 15;   end
if ~exist('rhoa', 'var'), rhoa  = 0.9;  end
if ~exist('rmi',  'var'), rmi   = 0.5;  end
if ~exist('nx',   'var'), nx    = 2;    end
if ~exist('ny',   'var'), ny    = 3;    end
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

u = randperm(n);
x = u(1:nx);
y = u(nx+1:nx+ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1 = ss_to_mvgc(A,C,K,V,x,y);

F2 = cpsd_to_mvgc(S,x,y,tol,maxi);

RE = 2*abs(F1-F2)/(abs(F1)+abs(F2));

fprintf('F1 = %8.6f\n',F1);
fprintf('F2 = %8.6f\n',F2);
fprintf('RE = %e\n\n',RE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = ss_to_smvgc(A,C,K,V,x,y,fres);

f2 = cpsd_to_smvgc(S,x,y,tol,maxi);

re = 2*maxabs(f1(:)-f2(:))/(maxabs(f1(:))+maxabs(f2(:)));

fprintf('re = %e\n\n',re);

om = sfreqs(fres);
gp_qplot(om,log([f1;f2])',{'SS','CPSD'},'set key outside\nset xlabel "angular frequency (rad)"\nset ylabel "spectral power (dB)" rot\nset xr[0:pi]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('n',         'var'), n         = 7;        end % number of variables
if ~exist('p',         'var'), p         = 4;        end % AR model order
if ~exist('rho',       'var'), rho       = 0.9;      end % spectral norm
if ~exist('w',         'var'), w         = 1.0;      end % AR decay weighting parameter (empty for no weighting)
if ~exist('g',         'var'), g         = 1.0;      end % residuals multi-information (0 for zero residuals correlation, empty for uniform random)
if ~exist('N',         'var'), N         = 1;        end % number of trials
if ~exist('m',         'var'), m         = 10000;    end % time series length
if ~exist('fs',        'var'), fs        = 200;      end % sample rate
if ~exist('fres',      'var'), fres      = [];       end % frequency resolution
if ~exist('window',    'var'), window    = [];       end % Welch window
if ~exist('noverlap',  'var'), noverlap  = [];       end % Welch overlap
if ~exist('tol',       'var'), tol       = [];       end % tolerance for Wilson's algorithm
if ~exist('maxiter',   'var'), maxiter   = [];       end % maximum iterations for Wilson's algorithm
if ~exist('mseed',     'var'), mseed     = 340219;   end % model random seed (0 to use current rng state)
if ~exist('dseed',     'var'), dseed     = 801123;   end % data random seed (0 to use current rng state)
if ~exist('gpterm',    'var'), gpterm    = 'x11';    end % Gnuplot terminal type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate VAR model

rstate = rng_seed(mseed);
A = var_rand(n,p,rho,w);
V = corr_rand(n,g);
rng_restore(rstate);

% generate time series

rstate = rng_seed(dseed);
X = var_to_tsdata(A,V,m,N);
rng_restore(rstate);

% estimate spectrum

fprintf('\nSpectral estimation ... ');tic;
[Se,f,nwobs,noobs,nwins] = tsdata_to_cpsd(X,fs,window,noverlap,fres);
et=toc;fprintf('%g seconds\n',et);
Pe = cpsd2power(Se);

fres = length(f)-1;

fprintf('\n');
fprintf('frequency resolution = %d\n',fres);
fprintf('number of windows    = %d\n',nwins);
fprintf('window  observations = %d\n',nwobs);
fprintf('overlap observations = %d\n',noobs);
fprintf('\n');

% calculate actual spectrum

[S,H] = var_to_cpsd(A,V,fres);
P = cpsd2power(S);

% plot estimated vs. actual spectrum

gplot_power({P,Pe},f,true,{'actual','estimated'},gpterm);

fprintf('\nrel. error: actual vs estimated spectrum = %g\n',relerr(S,Se));

% factorise actual spectrum

fprintf('\nWilson''s algorithm (actual spectrum) ... ');tic;
[Hf,Vf] = cpsd_specfac(S,tol,maxiter);
et=toc;fprintf('%g seconds\n',et);

Sf = var2cpsd(Hf,Vf);
fprintf('\nrel. error: actual vs factored actual spectrum = %g\n',relerr(S,Sf));

% factorise estimated spectrum

fprintf('\nWilson''s algorithm (estimated spectrum) ... ');tic;
[Hef,Vef] = cpsd_specfac(Se,tol,maxiter);
et=toc;fprintf('%g seconds\n',et);

Sef = var2cpsd(Hef,Vef);
fprintf('\nrel. error: estimated vs factored estimated spectrum = %g\n',  relerr(Se,Sef));
fprintf(  'rel. error: actual    vs factored estimated spectrum = %g\n\n',relerr(S,Sef));

function S = var2cpsd(H,V)

	[n,~,h] = size(H);
	S = zeros(n,n,h);
	L = chol(V,'lower');
	for k = 1:h
		M = H(:,:,k)*L;
		S(:,:,k) = M*M';
	end

end

function e = relerr(A,B)

	D = abs(B-A);
	M = abs(A);
	M(M <= 2*eps)=1; % Minimum detectable difference between x and a value close to x is O(x)*eps.
	E = D./M;
	e = mean(E(:));
%	e = max(E(:));

end

function gplot_power(P,f,dB,leg,gpterm)

	r = 0;
	if iscell(P)
		r = length(P);
		if isempty(leg), leg = cellstr(num2str((1:r)','var %d')); end
		for u = 1:r
			[h,n] = size(P{u}); % assume all the same
			if dB, P{u} = 20*log10(P{u}); end
			data{u} = [f P{u}];
		end
	else
		leg = 'var 1';
		[h,n] = size(P);
		if dB, P = 20*log10(P); end
		data = [f P];
	end

	gpstem = fullfile(tempdir,'gplot_power');

	gp_write(gpstem,data);

	gp = gp_open(gpstem,gpterm,[Inf 0.75]);
	fprintf(gp,'datfile = "%s.dat"\n',gpstem);
	fprintf(gp,'set xr[0:%g]\n',f(end));
	fprintf(gp,'set key outside\n');
	fprintf(gp,'set xlabel "Hz"\n');
	if dB
		fprintf(gp,'set ylabel "dB" norot\n');
	else
		fprintf(gp,'set ylabel "spectral density" rot\n');
	end
	fprintf(gp,'set multiplot layout %d,1\n',n);
	if r > 0
		for i = 1:n
			fprintf(gp,'plot \\\n');
			for u = 1:r
				fprintf(gp,'datfile i %d u 1:%d with lines ls %d t "%s", \\\n',u-1,1+i,u,leg{u});
			end
			fprintf(gp,'NaN not\n');
		end
	else
		for i = 1:n
			fprintf(gp,'plot datfile u 1:%d with lines ls 1 t "%s%\n',1+i,leg);
		end
	end
	fprintf(gp,'unset multiplot\n');
	gp_close(gp,gpstem,gpterm);

end

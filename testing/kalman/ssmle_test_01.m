%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('n',   'var'), n   = 5;    end
if ~exist('r',   'var'), r   = 7;    end
if ~exist('rho', 'var'), rho = 0.9;  end
if ~exist('g',   'var'), g   = 0.4;  end
if ~exist('m',   'var'), m   = 1000; end
if ~exist('vmm', 'var'), vmm = 30;   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verb = true;

[A,C,K] = iss_rand(n,r,rho);

V = corr_rand(n,g);

y = ss_to_tsdata(A,C,K,V,m);

varmo = tsdata_to_varmo(y,vmm,'LWR',[],true,true,false,'');

pf = 2*varmo

[moaic,mobic,mohqc,molrt] = tsdata_to_ssvarmo(y,pf,[],[],true,verb,'')

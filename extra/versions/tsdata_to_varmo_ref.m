function [moaic,mobic,mohqc,molrt] = tsdata_to_varmo_ref(X,q,regmode,alpha,motstat,pacf,plotit,verb)

if nargin < 3 || isempty(regmode), regmode = 'LWR';  end
if nargin < 4 || isempty(alpha),   alpha   = 0.005;  end
if nargin < 5 || isempty(motstat)
    ftest = true;
else
    switch lower(motstat)
        case 'f',    ftest = true;   % F-test form
        case 'chi2', ftest = false;  % chi2 test form
        otherwise,   error('unknown distribution (must be ''chi2'' or ''F'')');
    end
end

if nargin < 6 || isempty(pacf),    pacf    = 2;    end
if nargin < 7 || isempty(plotit),  plotit  = true; end
if nargin < 8 || isempty(verb),    verb    = 1;    end

if ~plotit, pacf = 0; end % no point if we're not going to display

n = size(X,1);

% get log-likelihoods at VAR model orders 1 .. q

[LL,Lk,Lm] = tsdata_to_mle_ref(X,q,regmode,verb>1); % log-likelihood, #free parameters. effective #observations

% calculate information criteria

[aic,bic,hqc] = infocrit(LL,Lk,Lm); % Akaike, Schwarz' Bayesian, Hannan-Quinn

% calculate optimal model orders according to information criteria (note: NaNs are ignored)

morder = (0:q)';
[~,idx] = min(aic); moaic = morder(idx);
[~,idx] = min(bic); mobic = morder(idx);
[~,idx] = min(hqc); mohqc = morder(idx);

% sequential likelihood ratio tests for optimal model order (Lutkephol)

lambda = 2*diff(LL); % LR test statistics
df = n*n;            % degrees of freedom (number of free parameters)
if ftest
    lrpval = 1-fcdf(lambda/df,df,Lm(2:end)-n*(1:q)'-1);
   %lrcrit = df*finv(1-lralpha,df,Lm(2:end)-n*(1:q)'-1);
else
    lrpval = 1-chi2cdf(lambda,df);
   %lrcrit = chi2inv(1-lralpha,df)*ones(q,1);
end
lralpha = alpha/q;   % Bonferroni correction on significance levels
hit = false;
for k = q:-1:1
    if lrpval(k) < lralpha, hit = true; break; end
end
if hit, molrt = k; else, molrt = 0; end

% partial autocorrelation

if pacf == 0
    R = [];
    Rcrit = NaN;
else
    if verb, fprintf('calculating partial autocorrelation... '); end
    if pacf == 1, fastpacf = false; elseif pacf == 2, fastpacf = true; else, error('bad ''pacf'' flag must be zero, 1 or 2'); end
    [R,Rm] = tsdata_to_pacf(X,q,fastpacf);  % partical autocorrelation function
    ccalpha = alpha/(n*n*q);                % Bonferroni correction on significance levels
    Rcrit = norminv(1-ccalpha/2)./sqrt(Rm); % 2-sided test
end

if verb > 1
    fprintf('\nBest model orders\n');
    fprintf('-----------------\n\n');
    fprintf('AIC : %2d',moaic); if moaic == q, fprintf(' *'); end; fprintf('\n');
    fprintf('BIC : %2d',mobic); if mobic == q, fprintf(' *'); end; fprintf('\n');
    fprintf('HQC : %2d',mohqc); if mohqc == q, fprintf(' *'); end; fprintf('\n');
    fprintf('LRT : %2d',molrt); if molrt == q, fprintf(' *'); end; fprintf('\n\n');
end

if plotit

    mo = morder(2:end);
    xo = 0.25;
    xlims = [1-xo q+xo];

    if moaic == q, waic = '*'; else, waic = ''; end
    if mobic == q, wbic = '*'; else, wbic = ''; end
    if mohqc == q, whqc = '*'; else, whqc = ''; end
    if molrt == q, wlrt = '*'; else, wlrt = ''; end

    aic1 = aic(2:end);
    bic1 = bic(2:end);
    hqc1 = hqc(2:end);
    saic = 0.05+0.95*(aic1-min(aic1))/(max(aic1)-min(aic1));
    sbic = 0.05+0.95*(bic1-min(bic1))/(max(bic1)-min(bic1));
    shqc = 0.05+0.95*(hqc1-min(hqc1))/(max(hqc1)-min(hqc1));
    if pacf, subplot(3,1,1); else, subplot(2,1,1); end
    plot(mo,[saic sbic shqc],'o-');
    grid on
    title('Information criteria');
    ylabel('information criterion (scaled)');
    xlabel('lags');
    legend(sprintf('AIC (opt = %d%c)',moaic,waic),sprintf('BIC (opt = %d%c)',mobic,wbic),sprintf('HQC (opt = %d%c)',mohqc,whqc));
    xlim(xlims);
    ylim([0 1.05]);
    %set(gca,'YTick',[]);

    lmin = eps; lpval = lrpval; lpval(lpval<=0) = lmin;
    if pacf, subplot(3,1,2); else, subplot(2,1,2); end
    semilogy(mo,lpval,'o-',mo,lralpha*ones(q,1),'r');
    grid on
    if ftest
        title(sprintf('Likelihood ratio F-test (\\alpha = %g)',alpha));
    else
        title(sprintf('Likelihood ratio \\chi^2-test (\\alpha = %g)',alpha));
    end
    ylabel('p-value');
    xlabel('lags');
    legend(sprintf('LRT (opt = %d%c)',molrt,wlrt),'location','southeast');
    xlim(xlims);
    ylim([lmin 100]);

    if pacf
        subplot(3,1,3);
        hold on
        for i = 1:n
            for j = 1:n
                plot(squeeze(R(i,j,:)),'o');
            end
        end
        plot(mo,[Rcrit -Rcrit],'r');
        hold off
        grid on
        title(sprintf('Partial autocorrelation (\\alpha = %g)',alpha));
        ylabel('correlation coefficient');
        xlabel('lags');
        xlim(xlims);
        ylim([-1 1]);
    end

	axes('Units','Normal');
	h = title(sprintf('VAR model order selection\n\n'),'FontSize',13);
	set(gca,'visible','off')
	set(h,'visible','on')

end

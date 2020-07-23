function X = simlfp(ROI,RCON,dt,simt,sett,usemex)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROI            struct (cell) array; ROI{r} is a struct with fields:
%
%    dect        vector, length = number of nodes, of decay times for nodes in R{r} (MILLISECONDS)
%    con         3-column matrix, number of rows = number of intra-ROI connections, specifying
%                a table of intra-ROI connections for ROI{r}: column 1 is the "to" node number,
%                column 2 the "from" node number, and column 3 the connection lag time (MILLISECONDS)
%
% RCON           5-column matrix, number of rows = number of inter-ROI connections, specifying
%                a table of inter-ROI connections: column 1 is the "to" ROI number, column 2
%                the "to" node number in the "to" ROI, column 3 is the "from" ROI number, column 4
%                the "from" node number in the "from" ROI, and  column 5 the connection lag time (MILLISECONDS)
%
% dt             Ornstein-Uhlenbeck time-slice increment (MILLISECONDS)
%
% simt           simulation time (SECONDS)
% sett           "settle" time, truncated at beginning (SECONDS)
%
% usemex         % use C (mex) version -- much faster
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrois = length(ROI);

for r = 1:nrois
	nnodes(r) = length(dect);
	lagtable



% convert millisecond times to seconds

dt    = dtm/1000;
lagt  = lagtm/1000;
mdect = mdectm/1000;
sdect = sdectm/1000;

% generate Gamma-distributed decay times

ga = (mdect^2)/(sdect^2); % Gamma shape parameter
gb = (sdect^2)/mdect;     % Gamma scale parameter
dect = gamrnd(ga,gb,n,1); % Gamma-distributed decay-times


X = oulag(1./dect,A,V,dt,lagt,simt,sett,usemex);

% Physical parameters
eps     = 0.000;
sigma   = 1.414222;
eta     = eps * sigma;

% Initial condition
pProp   = 1.190e-5;
n0      = 1.000;
r0      = (1 - pProp) * n0;
p0      = pProp * n0;

IC      = [log(n0); p0];

% Base switch rates
umin    = 1.200e-6;
v       = 0.100;
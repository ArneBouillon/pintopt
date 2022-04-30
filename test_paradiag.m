%% Tracking case
% Good example of eigenvalue diff between alpha=1 and alpha=-1:
%  d = 10, N = 100, Tend = 1e-4, gamma = 1e-4
d = 1;
N = 100;

Tend = 1e-4;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 10^-4;

K = diag(ones(d, 1)*2) - diag(ones(d-1, 1), 1) - diag(ones(d-1, 1), -1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);

yd = @(t) y0;
obj = Obj(ObjType.Tracking, gamma, yd);

K=112;

precinfo = Prec(PrecType.Tracking, struct('alpha', 1, 'test', true), obj);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, precinfo); toc

%% Terminal-cost case
d = 20;
N = 20;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d, 1)*2) - diag(ones(d-1, 1), 1) - diag(ones(d-1, 1), -1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCost, gamma, yT);

precinfo = Prec(PrecType.TerminalCost1, struct('alpha', -1, 'test', true), obj);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, precinfo); toc

%% Modified terminal-cost case
d = 10;
N = 10;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-2;

K = diag(ones(d, 1)*2) - diag(ones(d-1, 1), 1) - diag(ones(d-1, 1), -1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCostMod, gamma, yT);

precinfo = Prec(PrecType.TerminalCostMod, struct('alpha', -1, 'test', true), obj);
tic, [A,Ap] = paradiag(K, N, Tend, y0, obj, precinfo); toc

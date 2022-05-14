%% Tracking case
% Good example of eigenvalue diff between alpha=1 and alpha=-1:
%  d = 10, N = 100, Tend = 1e-4, gamma = 1e-4
d = 100;
N = 1000;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 10^-4;

K = diag(ones(d, 1)*2) - diag(ones(d-1, 1), 1) - diag(ones(d-1, 1), -1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;
K = sparse(K);

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);

yd = @(t) y0;
obj = Obj(ObjType.Tracking, gamma, yd);

precinfo = Prec(PrecType.Tracking, struct('alpha', -1, 'test', false), obj);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, precinfo); toc
return

%% Terminal-cost case
d = 1000;
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
K = sparse(K);

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCost, gamma, yT);

precinfo = Prec(PrecType.TerminalCost, struct('alpha', .0001, 'test', false), obj);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, precinfo); toc

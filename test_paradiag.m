%% Tracking case
d = 20;
N = 20;

Tend = .01;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 10^-3;

K = diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1);
K(1,end) = 1;
K(end,1) = 1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);

yd = @(t) y0;
obj = Obj(ObjType.Tracking, gamma, yd);

prectype = Prec(1, false, false);
[Y,L] = paradiag(K, N, Tend, y0, obj, prectype);

%% Terminal-cost case
d = 20;
N = 20;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1);
K(1,end) = 1;
K(end,1) = 1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCost, gamma, yT);

prectype = Prec(1, false, false);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, prectype); toc

%% Modified terminal-cost case
d = 20;
N = 20;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1);
K(1,end) = 1;
K(end,1) = 1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCostMod, gamma, yT);

prectype = Prec(1i, false, true);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, prectype); toc

%% Tracking case
d = 5;
N = 50;

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

[Y,L] = paradiag(K, N*10, Tend, y0, obj);

%% Terminal-cost case
d = 5;
N = 10;

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
yT = 0*yT;
obj = Obj(ObjType.TerminalCost, gamma, yT);

[Y,L] = paradiag(K, N*100, Tend, y0, obj);

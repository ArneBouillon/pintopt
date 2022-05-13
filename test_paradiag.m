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

gamma = 1e-4;

K = diag(ones(d, 1)*2) - diag(ones(d-1, 1), 1) - diag(ones(d-1, 1), -1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCostMod, gamma, yT);

precinfo = Prec(PrecType.TerminalCostMod, struct('alpha', -1, 'test', true), obj);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, precinfo); toc

%% Compare terminal-cost cases
d = 1;
N = 100;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d, 1)*2) - diag(ones(d-1, 1), 1) - diag(ones(d-1, 1), -1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;
K = 57;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCost, gamma, yT);
objmod = Obj(ObjType.TerminalCostMod, gamma, yT);
precinfo = Prec(PrecType.TerminalCost1, struct('alpha', -1, 'test', true), obj);
precinfomod = Prec(PrecType.TerminalCostMod, struct('alpha', -1, 'test', true), objmod);

[Y,L] = paradiag(K, N, Tend, y0, obj, precinfo);
[Ymod,Lmod] = paradiag(K, N, Tend, y0, objmod, precinfomod);

[Yacc,Lacc] = paradiag(K, N*10, Tend, y0, obj, precinfo);
Yacc=Yacc(:,1:10:end);
[Yaccmod,Laccmod] = paradiag(K, N*10, Tend, y0, objmod, precinfomod);
Yaccmod=Yaccmod(:,1:10:end);

f = @(t,yl) [
    -K*yl(1:d) - yl(d+1:end)/obj.gamma;
    K*yl(d+1:end);
];
bounds = @(yla, ylb) [
    yla(1:d) - y0;
    ylb(d+1:end) - ylb(1:d) + yT;
];

pts = N*50+1;
solinit = struct;
solinit.x = linspace(0, Tend, pts);
solinit.y = [y0*zeros(1, pts); y0*zeros(1, pts)];
sol = bvp5c(f, bounds, solinit, bvpset("RelTol", 1e-10));
ybvp=sol.y(1:d,1:50:end);

DT = Tend/N;
sh = DT*K;
gh = DT/gamma;

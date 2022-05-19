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

precinfo = Prec(PrecType.TerminalCost, struct('alpha', .0001, 'test', true), obj);
tic, [Y,L] = paradiag(K, N, Tend, y0, obj, precinfo); toc

%% Test of self-adjoint and non-self-adjoint equation against bvp5c
d = 10;
N = 100;

Tend = 1;
Ksymm = randn(d); Ksymm = Ksymm + Ksymm';
Knonsymm = randn(d);

gamma = 1e-0;

y0 = randn(d,1); yT = randn(d,1);
yd = @(t) t/Tend * yT + (Tend-t)/Tend * y0;

obj_tr = Obj(ObjType.Tracking, gamma, yd);
obj_tc = Obj(ObjType.TerminalCost, gamma, yT);
precinfo_tr = Prec(PrecType.Tracking, struct('alpha', -1, 'test', true), obj_tr);
precinfo_tc = Prec(PrecType.TerminalCost, struct('alpha', -1, 'test', true), obj_tc);

%% -> Symmetric tracking
f = @(t,yl) [
    -Ksymm*yl(1:d) - yl(d+1:end)/sqrt(gamma);
    (yd(t)-yl(1:d))/sqrt(gamma) + Ksymm'*yl(d+1:end);
];
bounds = @(yla, ylb) [
    yla(1:d) - y0;
    ylb(d+1:end);
];
solinit = struct;
solinit.x = linspace(0, Tend, N+1);
solinit.y = [y0*ones(1, N+1); y0*ones(1, N+1)];
sol_tr_symm = bvp5c(f, bounds, solinit, bvpset("RelTol", 1e-8));
[Y_tr_symm, L_tr_symm] = paradiag(Ksymm, N, Tend, y0, obj_tr, precinfo_tr, true);
plot(sol_tr_symm.x,sol_tr_symm.y(1:d,:)')
hold on, set(gca,'ColorOrderIndex',1)
plot(linspace(0,Tend,N+1), real(Y_tr_symm'), '.')

%% -> Non-symmetric tracking
f = @(t,yl) [
    -Knonsymm*yl(1:d) - yl(d+1:end)/sqrt(gamma);
    (yd(t)-yl(1:d))/sqrt(gamma) + Knonsymm'*yl(d+1:end);
];
bounds = @(yla, ylb) [
    yla(1:d) - y0;
    ylb(d+1:end);
];
solinit = struct;
solinit.x = linspace(0, Tend, N+1);
solinit.y = [y0*ones(1, N+1); y0*ones(1, N+1)];
sol_tr_nonsymm = bvp5c(f, bounds, solinit);
[Y_tr_nonsymm, L_tr_snonymm] = paradiag(Knonsymm, N, Tend, y0, obj_tr, precinfo_tr, true);
plot(sol_tr_nonsymm.x,sol_tr_nonsymm.y(1:d,:)')
hold on, set(gca,'ColorOrderIndex',1)
plot(linspace(0,Tend,N+1), real(Y_tr_nonsymm'), '.')

%% -> Symmetric terminal cost
f = @(t,yl) [
    -Ksymm*yl(1:d) - yl(d+1:end)/gamma;
    Ksymm'*yl(d+1:end);
];
bounds = @(yla, ylb) [
    yla(1:d) - y0;
    ylb(d+1:end)-ylb(1:d)+yT;
];
solinit = struct;
solinit.x = linspace(0, Tend, N+1);
solinit.y = [y0*ones(1, N+1); y0*ones(1, N+1)];
sol_tc_symm = bvp5c(f, bounds, solinit, bvpset("RelTol", 1e-8));
[Y_tc_symm, L_tc_symm] = paradiag(Ksymm, N, Tend, y0, obj_tc, precinfo_tc, true);
plot(sol_tc_symm.x,sol_tc_symm.y(1:d,:)')
hold on, set(gca,'ColorOrderIndex',1)
plot(linspace(0,Tend,N+1), real(Y_tc_symm'), '.')

%% -> Non-symmetric terminal cost
f = @(t,yl) [
    -Knonsymm*yl(1:d) - yl(d+1:end)/gamma;
    Knonsymm'*yl(d+1:end);
];
bounds = @(yla, ylb) [
    yla(1:d) - y0;
    ylb(d+1:end)-ylb(1:d)+yT;
];
solinit = struct;
solinit.x = linspace(0, Tend, N+1);
solinit.y = [y0*ones(1, N+1); y0*ones(1, N+1)];
sol_tc_nonsymm = bvp5c(f, bounds, solinit, bvpset("RelTol", 1e-8));
[Y_tc_nonsymm, L_tc_nonsymm] = paradiag(Knonsymm, N, Tend, y0, obj_tc, precinfo_tc, true);
plot(sol_tc_nonsymm.x,sol_tc_nonsymm.y(1:d,:)')
hold on, set(gca,'ColorOrderIndex',1)
plot(linspace(0,Tend,N+1), real(Y_tc_nonsymm'), '.')

%% Tracking case
d = 100;
N = 5;

Tend = .01;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);

yd = @(t) y0;
obj = Obj(ObjType.Tracking, gamma, yd);

prop_f = @(varargin) prop_ie_track(20, varargin{:});
prop_c = @(varargin) prop_ie_track(1, varargin{:});
% prop_f = @(varargin) prop_bvp5c_track(.000001, 50, varargin{:});
% prop_c = @(varargin) prop_bvp5c_track(.01, 5, varargin{:});
precinfo = Prec(PrecType.Tracking, struct('alpha', 1, 'test', false), obj);
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, Krylov.Specialized);

%% Terminal-cost case
d = 20;
N = 20;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCost, gamma, yT);

prop_f = @(varargin) prop_ie_tc(20, varargin{:});
prop_c = @(varargin) prop_ie_tc(1, varargin{:});
% prop_f = @(varargin) prop_bvp5c_tc(.000001, 50, varargin{:});
% prop_c = @(varargin) prop_bvp5c_tc(.001, 5, varargin{:});
precinfo = Prec(PrecType.TerminalCost2, struct('alpha', -1, 'test', true), obj);
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, Krylov.Generic);

%% Modified terminal-cost case
d = 10;
N = 10;

Tend = 1e-2;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);
yT = .5*exp(-100*(x-.25).^2)+.5*exp(-100*(x-.75).^2);
obj = Obj(ObjType.TerminalCostMod, gamma, yT);

prop_f = @(varargin) prop_ie_tc_mod(20, varargin{:});
prop_c = @(varargin) prop_ie_tc_mod(1, varargin{:});
% prop_f = @(varargin) prop_bvp5c_tc_mod(.000001, 500, varargin{:});
% prop_c = @(varargin) prop_bvp5c_tc_mod(.001, 5, varargin{:});
precinfo = Prec(PrecType.TerminalCostMod, struct('alpha', -1, 'test', false), obj);
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, Krylov.None);

%% Example for Krylov
d = 1000;
N = 5;

Tend = .5;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 1e-4;

K = diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1);
K(1,end) = -1;
K(end,1) = -1;
K = K/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);

yd = @(t) y0;
obj = Obj(ObjType.Tracking, gamma, yd);

DT = Tend / N;
gh = DT/sqrt(gamma);

prop_f = @(varargin) prop_ie_track(20, varargin{:});
prop_c = @(varargin) prop_ie_track(1, varargin{:});
precinfo = [];
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, Krylov.Specialized);

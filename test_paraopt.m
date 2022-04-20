%% Tracking case
d = 5;
N = 5;

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

prop_f = @(varargin) prop_ie_track(1000, varargin{:});
prop_c = @(varargin) prop_ie_track(1, varargin{:});
% prop_f = @(varargin) prop_bvp5c_track(.000001, 500, varargin{:});
% prop_c = @(varargin) prop_bvp5c_track(.001, 5, varargin{:});
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, [], Krylov.Generic);

%% Terminal-cost case
d = 50;
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
obj = Obj(ObjType.TerminalCost, gamma, yT);

prop_f = @(varargin) prop_ie_tc(100, varargin{:});
prop_c = @(varargin) prop_ie_tc(1, varargin{:});
% prop_f = @(varargin) prop_bvp5c_tc(.000001, 500, varargin{:});
% prop_c = @(varargin) prop_bvp5c_tc(.000001, 5, varargin{:});
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, [], Krylov.None);

%% Tracking case
d = 50;
N = 5;

Tend = .01;
xbegin = 0;
xend = 1;
dx = (xend - xbegin)/(d-1);

gamma = 10^-3;

A = -(diag(ones(d,1)*2)-diag(ones(d-1,1),1)-diag(ones(d-1,1),-1));
A(1,end) = 1;
A(end,1) = 1;
A = A/dx^2;

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);

yd = @(t) y0;
obj = Obj(ObjType.Tracking, gamma, yd);

prop_f = @(varargin) prop_ie_track(10, varargin{:});
prop_c = @(varargin) prop_ie_track(1, varargin{:});
[Y,L] = paraopt(A, N, Tend, y0, prop_f, prop_c, obj);

%% Terminal-cost case


%% Tracking case
d = 100;
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
K = sparse(K);

x = linspace(xbegin,xend,d)';
y0 = exp(-100*(x-.5).^2);

yd = @(t) y0;
obj = Obj(ObjType.Tracking, gamma, yd);

DT = Tend/N;
gh=DT/sqrt(gamma);
sh=DT*eig(K);

prop_f = @(varargin) prop_ie_track(10, varargin{:});
prop_c = @(varargin) prop_ie_track(1, varargin{:});
% prop_c = @prop_trap1_track;
% prop_f = @(varargin) prop_bvp5c_track(varargin{:}, .000001, 50);
% prop_c = @prop_bvp5c_track;
mp_c = MP_Track_IE1(K, obj, Tend/N);
precinfo = Prec(PrecType.Square, struct('alpha', -1, 'test', false), obj);
tic, [Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, SubEnh.Specialized, mp_c, 1e-6, 1e-4); toc

return
%% Terminal-cost case
d = 100;
N = 5;

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
% prop_f = @(varargin) prop_bvp5c_tc(varargin{:}, .000001, 50);
% prop_c = @prop_bvp5c_tc;
precinfo = Prec(PrecType.Square, struct('alpha', -1, 'test', false), obj);
mp_c = MP_TC_IE1(K, obj, Tend/N);
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, SubEnh.Specialized, mp_c);

%% Example for SubEnh
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
[Y,L] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, SubEnh.Specialized);

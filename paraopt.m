%% ParaOpt
% Perform ParaOpt to calculate the optimal control u(t) of the equation
%   (1)  y'(t) = Ay(t) + u(t),    y(0) = y0
% under a certain objective function.
%
% Parameters:
%  - A:         Matrix in the ODE (1) to control
%  - N:         Number of time intervals in which to split [0, Tend]
%  - Tend:      The end of the time interval [0, Tend] in which to solve
%               the problem
%  - y0:        The initial value in the ODE (1)
%  - prop_f:    The fine propagator to use
%  - prop_c:    The coarse propagator to use
%  - obj:       Info about the objective function (see Obj class)
%  - prec:      Info about the preconditioner to use
%                 Default: []
%  - krylov:    Information about whether to use Krylov-enhanced versions
%               of ParaOpt (see Krylov class)
%                 Default: Krylov.None
%  - tol:       Absolute tolerance of the ParaOpt iteration
%                 Default: 10^-8
%  - Y0:        Initial guesses for the Y variables
%                 Default: random
%  - L0:        Initial guesses for the L variables
%                 Default: random
%
% Return values:
%  - Y:      Discretised solution to (1)
%  - L:      Discretised adjoint
%  - k:      Number of (outer) ParaOpt iterations
%
% Internal data layout:
%  - Size of systems: - 2*d*(N-1) for tracking-type objectives
%                     - 2*d*N     for terminal-cost objectives
%  - Y :: (d,N+1):    Y_0, Y_1, ..., Y_N at second indices 1, 2, ..., N+1
%  - L :: (d,N+1):    L_0, L_1, ..., L_N at second indices 1, 2, ..., N+1
%
function [Y,L,k] = paraopt(A, N, Tend, y0, prop_f, prop_c, obj, prec, ...
                           krylov, tol, Y0, L0)
    d = size(A,1);
    DT = Tend / N;

    rng(1337)
    if ~exist('prec', 'var') || isempty([]), prec = []; end
    if ~exist('krylov', 'var') || isempty(krylov), krylov = Krylov.None; end
    if ~exist('tol', 'var') || isempty(tol), tol = 10^-8; end
    if ~exist('Y0', 'var') || isempty(Y0), Y0 = randn(d, N+1); end
    if ~exist('L0', 'var') || isempty(L0), L0 = randn(d, N+1); end

    [Y,L] = init_YL(Y0, L0, obj, y0);

    P00 = []; Q00 = []; S = []; SP = []; SQ = [];
    if krylov.any
        P00 = zeros(d, N); Q00 = zeros(d, N);
        for n=1:N
            [P,Q] = prop_f(zeros(d,1), zeros(d,1), (n-1)*DT, n*DT, obj, A, false);
            P00(:,n) = P; Q00(:,n) = Q;
        end
    end

    k = 0; gmresiter = 0;
    while true
        Ps = zeros(d, N); Qs = zeros(d, N);
        for n=1:N
            [P,Q] = prop_f(Y(:,n), L(:,n+1), (n-1)*DT, n*DT, obj, A, false);
            Ps(:,n) = P; Qs(:,n) = Q;
        end
        [Y,L] = fill_in(Y, L, Ps, Qs, obj);

        switch krylov
            case Krylov.Generic
                [S, SP, SQ] = add_orth(S, SP, SQ, [Y(:,1:N); L(:,2:end)], Ps-P00, Qs-Q00);
            case Krylov.Specialized
                switch obj.type
                    case ObjType.Tracking
                        [S, SP, SQ] = add_orth(S, SP, SQ, [[Y(:,1:N); L(:,2:end)] [L(:,2:end); -Y(:,1:N)]], [Ps-P00 Qs-Q00], [Qs-Q00 -(Ps-P00)]);
                    case ObjType.TerminalCost
                        [S, SP, SQ] = add_orth(S, SP, SQ, [[Y(:,1:N); L(:,2:end)] [Y(:,1:N)+L(:,2:end); L(:,2:end)]], [Ps-P00 Ps+Qs-P00-Q00], [Qs-Q00 Qs-Q00]);
                end
            case Krylov.None
        end
        F = get_F(Y, L, Ps, Qs, obj);
        nrm = norm(F);
        disp(['Iteration ' num2str(k) ': ' num2str(nrm) ' in ' num2str(gmresiter(end)) ' GMRES iterations'])
        if nrm < tol, break, end
        k = k + 1;

        apply_jac_fun = get_apply_jac_fun(A, @krylov_prop_c, obj, N, DT);
        [delta,~,~,gmresiter] = gmres(apply_jac_fun, -F, [], [], numel(F), prec);

        dY = reshape(delta(1:numel(delta)/2), d, []);
        dL = reshape(delta(numel(delta)/2+1:end), d, []);
        Y(:,2:1+size(dY,2)) = Y(:,2:1+size(dY,2)) + dY;
        L(:,2:1+size(dL,2)) = L(:,2:1+size(dL,2)) + dL;
    end
    
    function [P,Q] = krylov_prop_c(dy, dl, varargin)
        if krylov.any
            comps = S' * [dy; dl];
            ylapprox = S * comps; yapprox = ylapprox(1:d); lapprox = ylapprox(d+1:end);
            [P,Q] = prop_c(dy - yapprox, dl - lapprox, varargin{:});
            P = P + SP * comps;
            Q = Q + SQ * comps;
        else
            [P,Q] = prop_c(dy, dl, varargin{:});
        end
    end
end

function [Y,L] = init_YL(Y0, L0, obj, y0)
    Y = Y0; L = L0;
    Y(:,1) = y0;
    L(:,1) = NaN;
    switch obj.type
        case ObjType.Tracking, Y(:,end) = NaN; L(:,end) = 0;
        case ObjType.TerminalCost
    end
end

function [Y,L] = fill_in(Y, L, Ps, Qs, obj)
    L(:,1) = Qs(:,1);
    switch obj.type
        case ObjType.Tracking, Y(:,end) = Ps(:,end);
        case ObjType.TerminalCost
    end
end

function F = get_F(Y, L, Ps, Qs, obj)
    switch obj.type
        case ObjType.Tracking
            F = [
                reshape(Y(:,2:end-1) - Ps(:,1:end-1), [], 1);
                reshape(L(:,2:end-1) - Qs(:,2:end), [], 1);
            ];
        case ObjType.TerminalCost
            F = [
                reshape(Y(:,2:end)   - Ps(:,1:end), [], 1);
                reshape(L(:,2:end-1) - Qs(:,2:end), [], 1);
                L(:,end) - Y(:,end) + obj.y_T;
            ];
    end
end

function apply_jac_fun = get_apply_jac_fun(A, prop_c, obj, N, DT)
    switch obj.type
        case ObjType.Tracking, apply_jac_fun = @(delta) apply_jac_track(delta, A, prop_c, obj, N, DT);
        case ObjType.TerminalCost, apply_jac_fun = @(delta) apply_jac_tc(delta, A, prop_c, obj, N, DT);
    end
end

function res = apply_jac_track(delta, A, prop_c, obj, N, DT)
    dY = reshape(delta(1:numel(delta)/2), [], N-1);
    dL = reshape(delta(numel(delta)/2+1:end), [], N-1);

    d = size(dY,1);

    dY0 = zeros(d,1);
    dLend = zeros(d,1);

    res = delta;
    for n=1:N
        if n == 1, dy = dY0; else, dy = dY(:,n-1); end
        if n == N, dl = dLend; else, dl = dL(:,n); end
        [P,Q] = prop_c(dy, dl, (n-1)*DT, n*DT, obj, A, true);
        if n < N
            res((n-1)*d+1:n*d) = res((n-1)*d+1:n*d) - P;
        end
        if n > 1
            res((N-1)*d+(n-2)*d+1:(N-1)*d+(n-1)*d) = res((N-1)*d+(n-2)*d+1:(N-1)*d+(n-1)*d) - Q;
        end
    end
end

function res = apply_jac_tc(delta, A, prop_c, obj, N, DT)
    dY = reshape(delta(1:numel(delta)/2), [], N);
    dL = reshape(delta(numel(delta)/2+1:end), [], N);

    d = size(dY,1);

    dY0 = zeros(d,1);

    res = delta;
    for n=1:N
        if n == 1, dy = dY0; else, dy = dY(:,n-1); end
        dl = dL(:,n);
        [P,Q] = prop_c(dy, dl, (n-1)*DT, n*DT, obj, A, true);
        res((n-1)*d+1:n*d) = res((n-1)*d+1:n*d) - P;
        if n > 1
            res(N*d+(n-2)*d+1:N*d+(n-1)*d) = res(N*d+(n-2)*d+1:N*d+(n-1)*d) - Q;
        end
    end

    res(end-d+1:end) = res(end-d+1:end) - dY(:,end);
end

function [S, SP, SQ] = add_orth(S, SP, SQ, newS, newSP, newSQ)
    for i=1:size(newS,2)
        s = newS(:,i); sp = newSP(:,i); sq = newSQ(:,i);

        for ii=1:2 % Re-orthogonalisation
            if S
                comps = S'*s;
                s = s - S*comps; sp = sp - SP*comps; sq = sq - SQ*comps;
            end

            scale = norm(s);
            if scale < sqrt(eps), break, end
            s = s / scale; sp = sp / scale; sq = sq / scale;

            if ii==2, S = [S s]; SP = [SP sp]; SQ = [SQ sq]; end
        end
    end
end

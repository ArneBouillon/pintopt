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
%  - compextra: Whether to compute the components in Y and L that aren't
%               needed to compute the rest, instead of setting them to NaN
%                 Default: true
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
                           krylov, tol, compextra, Y0, L0)
    d = size(A,1);

    if ~exist('prec', 'var'), prec = []; end
    if ~exist('krylov', 'var'), krylov = false; end
    if ~exist('tol', 'var'), tol = 10^-8; end
    if ~exist('compextra', 'var'), compextra = true; end
    if ~exist('Y0', 'var'), Y0 = randn(d, N+1); end
    if ~exist('L0', 'var'), L0 = randn(d, N+1); end

    [Y,L] = init_YL(Y0, L0, obj, y0);
    
    if krylov.any
        % TODO: Calculate 00
    end

    k = 0;
    nrm = +inf;
    while nrm > tol
        k = k + 1;

        Ps = zeros(d, N); Qs = zeros(d, N);
        for n=1:N
            [P,Q] = prop_f(A, Y(:,n), L(:,n+1));
            Ps(:,n) = P; Qs(:,n) = Q;
        end

        F = get_F(Y, L, Ps, Qs, obj);
        apply_jac_fun = get_apply_jac_fun(A, prop_c, obj, krylov, N);
        
        [delta,~,~,iter] = gmres(apply_jac_fun, -F, [], [], numel(F), prec);
        nrm = norm(F);
        
        disp(['Iteration ' num2str(k) ': ' num2str(nrm) ' in ' num2str(iter(end)) ' GMRES iterations'])

        dY = reshape(delta(1:numel(delta)/2), d, []);
        dL = reshape(delta(numel(delta)/2+1:end), d, []);
        Y(:,2:1+size(dY,2)) = Y(:,2:1+size(dY,2)) + dY;
        L(:,2:1+size(dL,2)) = L(:,2:1+size(dL,2)) + dL;
    end
    
    if compextra
        % TODO: Comp extra
    end
end

function [Y,L] = init_YL(Y0, L0, obj, y0)
    Y = Y0; L = L0;
    Y(:,1) = y0;
    L(:,1) = NaN;
    switch obj.type
        case ObjType.Tracking, Y(:,end) = NaN; L(:,end) = NaN;
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
                L(:,end) - Y(:,end) + obj.y_d;
            ];
    end
end

function apply_jac_fun = get_apply_jac_fun(A, prop_c, obj, krylov, N)
    switch obj.type
        case ObjType.Tracking, apply_jac_fun = @(delta) apply_jac_track(delta, A, prop_c, krylov, N);
        case ObjType.TerminalCost, apply_jac_fun = @(delta) apply_jac_tc(delta, A, prop_c, krylov, N);
    end
end

function res = apply_jac_track(delta, A, prop_c, krylov, N)
    dY = reshape(delta(1:numel(delta)/2), [], N-1);
    dL = reshape(delta(numel(delta)/2+1:end), [], N-1);

    d = size(dY,1);

    dY0 = zeros(d,1);
    dLend = zeros(d,1);

    switch krylov
        % TODO
    end

    res = delta;

    for n=1:N
        if n == 1
            dy = dY0;
        else
            dy = dY(:,n-1);
        end
        [P,Q] = prop_c(A, dy, zeros(d,1));
        if n < N-1
            res(n*d+1:(n+1)*d) = res(n*d+1:(n+1)*d) - P;
        end
        if n > 0
            res((N-1)*d+(n-1)*d+1:(N-1)*d+n*d) = res((N-1)*d+(n-1)*d+1:(N-1)*d+n*d) - Q;
        end
    end
    for n=1:N
        if n == N
            dl = dLend;
        else
            dl = dL(:,n);
        end
        [P,Q] = prop_c(A, zeros(d,1), dl);
        if n < N
            res((n-1)*d+1:n*d) = res((n-1)*d+1:n*d) - P;
        end
        if n > 1
            res((N-1)*d+(n-2)*d+1:(N-1)*d+(n-1)*d) = res((N-1)*d+(n-2)*d+1:(N-1)*d+(n-1)*d) - Q;
        end
    end
end

function apply_jac_tc(delta, A, prop_c, krylov, N)
    % TODO
    raise
end

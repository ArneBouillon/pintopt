%% ParaOpt
% Perform ParaOpt to calculate the optimal control u(t) of the equation
%   (1)  y'(t) = -Ky(t) + u(t),    y(0) = y0
% under a certain objective function.
%
% Parameters:
%  - K:         Matrix in the ODE (1) to control
%  - N:         Number of time intervals in which to split [0, Tend]
%  - Tend:      The end of the time interval [0, Tend] in which to solve
%               the problem
%  - y0:        The initial value in the ODE (1)
%  - prop_f:    The fine propagator to use
%  - prop_c:    The coarse propagator to use
%  - obj:       Info about the objective function (see Obj class)
%  - precinfo:  Info about the preconditioner to use
%                 Default: No preconditioner
%  - krylov:    Information about whether to use Krylov-enhanced versions
%               of ParaOpt (see Krylov class)
%                 Default: Krylov.None
%  - tol:       Absolute tolerance of the ParaOpt iteration
%                 Default: 10^-8
%  - gmrestol:  Relative residual tolerance of the inner GMRES solver
%                 Default: 10^-3
%  - Y0:        Initial guesses for the Y variables
%                 Default: random
%  - L0:        Initial guesses for the L variables
%                 Default: random
%
% Return values:
%  - Y:      Discretised solution to (1)
%  - L:      Discretised adjoint
%  - k:      Number of (outer) ParaOpt iterations
%  - res:    An array containing the residual norm in each iteration
%
% Internal data layout:
%  - Size of systems: - 2*d*(N-1) for tracking-type objectives
%                     - 2*d*N     for terminal-cost objectives
%  - Y :: (d,N+1):    Y_0, Y_1, ..., Y_N at second indices 1, 2, ..., N+1
%  - L :: (d,N+1):    L_0, L_1, ..., L_N at second indices 1, 2, ..., N+1
%
function [Y,L,k,res] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, ...
                               krylov, tol, gmrestol, Y0, L0)
    d = size(K,1);
    DT = Tend / N;

    rng(1337)
    if ~exist('precinfo', 'var') || isempty(precinfo), precinfo = []; end
    if ~exist('krylov', 'var') || isempty(krylov), krylov = Krylov.None; end
    if ~exist('tol', 'var') || isempty(tol), tol = 10^-8; end
    if ~exist('gmrestol', 'var') || isempty(gmrestol), gmrestol = 10^-3; end
    if ~exist('Y0', 'var') || isempty(Y0), Y0 = randn(d, N+1); end
    if ~exist('L0', 'var') || isempty(L0), L0 = randn(d, N+1); end

    [Y,L] = init_YL(Y0, L0, obj, y0);

    P00 = []; Q00 = []; S = []; SP = []; SQ = [];
    if krylov.any
        P00 = zeros(d, N); Q00 = zeros(d, N);
        for n=1:N
            [P,Q] = prop_f(zeros(d,1), zeros(d,1), (n-1)*DT, n*DT, obj, K, false);
            P00(:,n) = P; Q00(:,n) = Q;
        end
    end

    k = 0; gmresiter = 0; res = []; flag = NaN; relres = NaN;
    while true
        Ps = zeros(d, N); Qs = zeros(d, N);
        for n=1:N
            [P,Q] = prop_f(Y(:,n), L(:,n+1), (n-1)*DT, n*DT, obj, K, false);
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
                    case ObjType.TerminalCostMod
                        error 'Not implemented'
                end
            case Krylov.None
        end
        F = get_F(Y, L, Ps, Qs, obj);
        nrm = norm(F);
        res = [res nrm];
        disp(['Iteration ' num2str(k) ': ' num2str(nrm) ' in ' num2str(gmresiter(end)) ' GMRES iterations (flag=' num2str(flag) ', relres=' num2str(relres) ')'])
        if nrm < tol, break, end
        k = k + 1;

        apply_jac_fun = get_apply_jac_fun(K, @krylov_prop_c, obj, N, DT);
        prec = get_prec(apply_jac_fun, K, @krylov_prop_c, obj, N, DT, precinfo);

%         j = zeros(numel(F)); p = zeros(numel(F));
%         for i=1:numel(F)
%             disp([num2str(i) '/' num2str(numel(F))])
%             v = zeros(numel(F),1);
%             v(i) = 1;
%             j(:,i) = apply_jac_fun(v);
%             p(:,i) = prec(v);
%         end
%         n = size(j,1);
%         M = [speye(n/2,n/2) sparse(n/2,n/2);speye(n/2,n/2) speye(n/2,n/2)];
% %         figure, surf(real(p*j))
% %         figure, imshow(abs(j))
% %         figure, imshow(abs(inv(p)))
% %         figure, imshow(abs(inv(p)-j))
% %         figure, surf(real(inv(p)))
% %         figure, imshow(abs(inv(p)-prc))
% %         figure, surf(abs(p*j))
% %         j((1:15)+(N-1)*d,1:15),i=inv(p);i((1:15)+(N-1)*d,1:15)
% %         j((1:15)+(N-1)*d,1:15)./p((1:15)+(N-1)*d,1:15)
% % %         m=p*j;m(1:10,1:10)
% % %         norm(abs(p*j-eye(numel(F))))
%         condest(p*j)
%         Y=j;L=p;return
%         pause

        [delta,flag,relres,gmresiter] = gmres(apply_jac_fun, -F, [], gmrestol, 100, prec);

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
        case {ObjType.Tracking, ObjType.TerminalCostMod}, Y(:,end) = NaN; L(:,end) = 0;
        case ObjType.TerminalCost
    end
end

function [Y,L] = fill_in(Y, L, Ps, Qs, obj)
    L(:,1) = Qs(:,1);
    switch obj.type
        case {ObjType.Tracking, ObjType.TerminalCostMod}, Y(:,end) = Ps(:,end);
        case ObjType.TerminalCost
    end
end

function F = get_F(Y, L, Ps, Qs, obj)
    switch obj.type
        case {ObjType.Tracking, ObjType.TerminalCostMod}
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

function prec = get_prec(A, K, krylov_prop_c, obj, N, DT, precinfo)
    d = size(K, 1);

    if isempty(precinfo)
        prec = [];
        return
    end

    if precinfo.test
        prec = get_test_prec(A, d, N, obj, precinfo);
        return
    end

    prec = @(vec) prec_fun(vec, K, krylov_prop_c, obj, N, DT, precinfo);
end

function prec = get_test_prec(Afun, d, N, obj, precinfo)
    switch obj.type
        case ObjType.TerminalCost, n = d*N*2;
        case {ObjType.Tracking, ObjType.TerminalCostMod}, n = d*(N-1)*2;
    end

    A = sparse(n, n);
    for i=1:n
        v = zeros(n, 1);
        v(i) = 1;
        A(:,i) = Afun(v);
    end

    switch precinfo.type
        case {PrecType.Tracking, PrecType.TerminalCostMod}
            error 'Not implemented'
        case PrecType.TerminalCost1
            prec = A;
            prec(n/2+1:end,1:n/2) = 0;
            prec(1:d,n/2-d+1:n/2) = precinfo.alpha * prec(d+1:2*d,1:d);
            prec(n-d+1:n,n/2+1:n/2+d) = precinfo.alpha' * prec(n/2+1:n/2+d,n/2+d+1:n/2+2*d);
        case PrecType.TerminalCost2
            M = [speye(n/2,n/2) sparse(n/2,n/2);speye(n/2,n/2) speye(n/2,n/2)];
            prec = M*A*M;

            prec(1:d,n/2-2*d+1:n/2-d) = precinfo.alpha * prec(d+1:2*d,1:d);
            prec(end-2*d+1:end-d,n/2+1:n/2+d) = precinfo.alpha' * prec(n/2+1:n/2+d,n/2+d+1:n/2+2*d);
            prec(end-d+1:end,end-d+1:end) = speye(d);
            prec(end-d+1:end,1:n/2) = 0;

            prec(end-2*d+1:end-d,n/2-d+1:n/2) = 0;
            prec(end-2*d+1:end-d,1:d) = precinfo.alpha * prec(n/2+1:n/2+d,d+1:2*d);
            prec(n/2+1:n/2+d,n/2-2*d+1:n/2-d) = precinfo.alpha' * prec(n/2+d+1:n/2+2*d,1:d);

            prec = M\prec/M;
    end
end

function res = prec_fun(vec, K, krylov_prop_c, obj, N, DT, precinfo)
    d = size(K, 1);

    switch precinfo.type % Do any preliminary work
        case PrecType.TerminalCost2
            n = numel(vec);
            vec(n/2+1:end) = vec(n/2+1:end) + vec(1:n/2);

            lend = vec(end-d+1:end);
            [P,Q] = krylov_prop_c(zeros(d,1), lend, 0, DT, obj, K, true);
            vec(n/2-d+1:n/2) = vec(n/2-d+1:n/2) + P;
            vec(end-2*d+1:end-d) = vec(end-2*d+1:end-d) + Q;

            byend = vec(n/2-d+1:n/2);
            vec = [vec(1:n/2-d); vec(n/2+1:end-d)];
    end
    
    switch precinfo.type % Set the size of the system
        case PrecType.TerminalCost1, M = N;
        case {PrecType.Tracking, PrecType.TerminalCostMod, PrecType.TerminalCost2}, M = N - 1;
    end

    % Step 1: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * precinfo.alpha^((m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * 1/(precinfo.alpha^((m-1)/M))';
    end

    % Step 2: FF
    vec = reshape(vec, d, 2*M);
    vec = [ifft(vec(:,1:M).'); ifft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)*sqrt(M);

    % Step 3: Solve systems
    D = zeros(M,1); D(2) = precinfo.alpha^(1/M); D = M*ifft(D);
    D2 = zeros(M,1); D2(2) = precinfo.alpha^(1/M); D2(end) = precinfo.alpha^((M-1)/M); D2 = M*ifft(D2);
    for m=1:M
        [sol,~,~,~] = gmres(...
            @subsys, ...
            [vec((m-1)*d+1:m*d); vec(M*d+(m-1)*d+1:M*d+m*d)], ...
            [], 1e-10, 2*d...
        );
        vec((m-1)*d+1:m*d) = sol(1:d);
        vec(M*d+(m-1)*d+1:M*d+m*d) = sol(d+1:end);
    end

    % Step 4: F'
    vec = reshape(vec, d, 2*M);
    vec = [fft(vec(:,1:M).'); fft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)/sqrt(M);

    % Step 5: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * precinfo.alpha^(-(m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * (precinfo.alpha^((m-1)/M))';
    end

    switch precinfo.type % Do any concluding work
        case PrecType.TerminalCost2
            [P,~] = krylov_prop_c(vec(n/2-2*d+1:n/2-d), zeros(d,1), m*DT, (m+1)*DT, obj, K, true);
            [yend,~,~,~] = gmres(@subsys2, byend + P, [], [], d);
            vec = [vec(1:n/2-d); yend; vec(n/2-d+1:end); lend];
            vec(n/2+1:end) = vec(n/2+1:end) + vec(1:n/2);
    end

    res = vec;

    function prd = subsys(v)
        switch precinfo.type
            case PrecType.TerminalCost2
                prd = v;
                prd(d+1:end) = prd(d+1:end) + 2*v(1:d);
                [P,~] = krylov_prop_c(D(m)*v(1:d), (v(1:d)+v(d+1:end))*(1-precinfo.diag), m*DT, (m+1)*DT, obj, K, true);
                prd(1:d) = prd(1:d) - P;
                [P,~] = krylov_prop_c(D(m)'*v(d+1:end) + D2(m)*v(1:d), v(1:d)+v(d+1:end), m*DT, (m+1)*DT, obj, K, true);
                prd(d+1:end) = prd(1+d:end) - P;
            case {PrecType.Tracking, PrecType.TerminalCost1, PrecType.TerminalCostMod}
                prd = v;
                [P,~] = krylov_prop_c(D(m)*v(1:d), v(d+1:end)*(1-precinfo.diag), m*DT, (m+1)*DT, obj, K, true);
                prd(1:d) = prd(1:d) - P;
                [~,Q] = krylov_prop_c(v(1:d)*(1-precinfo.diag), D(m)'*v(d+1:end), m*DT, (m+1)*DT, obj, K, true);
                prd(d+1:end) = prd(d+1:end) - Q;
        end
    end

    function prd = subsys2(v)
        prd = v;
        [P,~] = krylov_prop_c(zeros(d, 1), -v, m*DT, (m+1)*DT, obj, K, true);
        prd = prd + P;
    end
end

function apply_jac_fun = get_apply_jac_fun(K, prop_c, obj, N, DT)
    switch obj.type
        case ObjType.Tracking, apply_jac_fun = @(delta) apply_jac_track(delta, K, prop_c, obj, N, DT);
        case ObjType.TerminalCost, apply_jac_fun = @(delta) apply_jac_tc(delta, K, prop_c, obj, N, DT);
        case ObjType.TerminalCostMod, apply_jac_fun = @(delta) apply_jac_tc_mod(delta, K, prop_c, obj, N, DT);
    end
end

function res = apply_jac_track(delta, K, prop_c, obj, N, DT)
    dY = reshape(delta(1:numel(delta)/2), [], N-1);
    dL = reshape(delta(numel(delta)/2+1:end), [], N-1);

    d = size(dY,1);

    dY0 = zeros(d,1);
    dLend = zeros(d,1);

    res = delta;
    for n=1:N
        if n == 1, dy = dY0; else, dy = dY(:,n-1); end
        if n == N, dl = dLend; else, dl = dL(:,n); end
        [P,Q] = prop_c(dy, dl, (n-1)*DT, n*DT, obj, K, true);
        if n < N
            res((n-1)*d+1:n*d) = res((n-1)*d+1:n*d) - P;
        end
        if n > 1
            res((N-1)*d+(n-2)*d+1:(N-1)*d+(n-1)*d) = res((N-1)*d+(n-2)*d+1:(N-1)*d+(n-1)*d) - Q;
        end
    end
end

function res = apply_jac_tc(delta, K, prop_c, obj, N, DT)
    dY = reshape(delta(1:numel(delta)/2), [], N);
    dL = reshape(delta(numel(delta)/2+1:end), [], N);

    d = size(dY,1);

    dY0 = zeros(d,1);

    res = delta;
    for n=1:N
        if n == 1, dy = dY0; else, dy = dY(:,n-1); end
        dl = dL(:,n);
        [P,Q] = prop_c(dy, dl, (n-1)*DT, n*DT, obj, K, true);
        res((n-1)*d+1:n*d) = res((n-1)*d+1:n*d) - P;
        if n > 1
            res(N*d+(n-2)*d+1:N*d+(n-1)*d) = res(N*d+(n-2)*d+1:N*d+(n-1)*d) - Q;
        end
    end

    res(end-d+1:end) = res(end-d+1:end) - dY(:,end);
end

function res = apply_jac_tc_mod(delta, K, prop_c, obj, N, DT)
    % These are the same
    res = apply_jac_track(delta, K, prop_c, obj, N, DT);
end

function [S, SP, SQ] = add_orth(S, SP, SQ, newS, newSP, newSQ)
    for i=1:size(newS,2)
        s = newS(:,i); sp = newSP(:,i); sq = newSQ(:,i);

        for ii=1:2 % Re-orthogonalisation
            if ~isempty(S)
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

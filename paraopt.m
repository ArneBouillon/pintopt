%% Linear ParaOpt
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
%  - subenh:    Information about whether to use subspace-enhanced versions
%               of ParaOpt (see SubEnh class)
%                 Default: SubEnh.None
%  - tol:       Absolute tolerance of the ParaOpt iteration
%                 Default: 10^-8
%  - gmrestol:  Relative residual tolerance of the inner GMRES solver
%                 Default: 10^-3
%  - gmresmxit: Maximum number of iterations of the inner GMRES solver
%                 Default: 50
%  - silent:    Mute all output to the console
%                 Default: false
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
%  - kgmres: An array containing the number of GMRES iterations in each
%            outer ParaOpt iteration
%
% Internal data layout:
%  - Size of systems: - 2*d*(N-1) for tracking-type objectives
%                     - 2*d*N     for terminal-cost objectives
%  - Y :: (d,N+1):    Y_0, Y_1, ..., Y_N at second indices 1, 2, ..., N+1
%  - L :: (d,N+1):    L_0, L_1, ..., L_N at second indices 1, 2, ..., N+1
%
function [Y,L,k,res,kgmres] = paraopt(K, N, Tend, y0, prop_f, prop_c, obj, precinfo, ...
                                      subenh, mp_c, tol, gmrestol, gmresmxit, ...
                                      silent, Y0, L0)
    d = size(K,1);
    DT = Tend / N;

    rng(1337)
    if ~exist('precinfo', 'var') || isempty(precinfo), precinfo = []; end
    if ~exist('subenh', 'var') || isempty(subenh), subenh = SubEnh.None; end
    if ~exist('mp_c', 'var') || isempty(mp_c), mp_c = []; end
    if ~exist('tol', 'var') || isempty(tol), tol = 10^-8; end
    if ~exist('gmrestol', 'var') || isempty(gmrestol), gmrestol = 10^-3; end
    if ~exist('gmresmxit', 'var') || isempty(gmresmxit), gmresmxit = 50; end
    if ~exist('silent', 'var') || isempty(silent), silent = false; end
    if ~exist('Y0', 'var') || isempty(Y0), Y0 = randn(d, N+1); end
    if ~exist('L0', 'var') || isempty(L0), L0 = randn(d, N+1); end

    [Y,L] = init_YL(Y0, L0, obj, y0);

    P00 = []; Q00 = []; Pc00 = []; Qc00 = []; S = []; SP = []; SQ = [];
    if subenh.any
        P00 = zeros(d, N); Q00 = zeros(d, N); Pc00 = zeros(d, N); Qc00 = zeros(d, N);
        for n=1:N
            [P,Q] = prop_f(zeros(d,1), zeros(d,1), (n-1)*DT, n*DT, obj, K, false);
            P00(:,n) = P; Q00(:,n) = Q;
            [P,Q] = prop_c(zeros(d,1), zeros(d,1), (n-1)*DT, n*DT, obj, K, false);
            Pc00(:,n) = P; Qc00(:,n) = Q;
        end
    end

    k = 0; gmresiter = 0; res = []; kgmres = []; flag = NaN; relres = NaN;
    while true
        Ps = zeros(d, N); Qs = zeros(d, N);
        for n=1:N
            [P,Q] = prop_f(Y(:,n), L(:,n+1), (n-1)*DT, n*DT, obj, K, false);
            Ps(:,n) = P; Qs(:,n) = Q;
        end
        [Y,L] = fill_in(Y, L, Ps, Qs, obj);

        switch subenh
            case SubEnh.Generic
                [S, SP, SQ] = add_orth(S, SP, SQ, [Y(:,1:N); L(:,2:end)], Ps-P00, Qs-Q00);
            case SubEnh.Specialized
                switch obj.type
                    case ObjType.Tracking
                        [S, SP, SQ] = add_orth(S, SP, SQ, [[Y(:,1:N); L(:,2:end)] [L(:,2:end); -Y(:,1:N)]], [Ps-P00 Qs-Q00], [Qs-Q00 -(Ps-P00)]);
                    case ObjType.TerminalCost
                        [S, SP, SQ] = add_orth(S, SP, SQ, [[Y(:,1:N); L(:,2:end)] [Y(:,1:N)+L(:,2:end); L(:,2:end)]], [Ps-P00 Ps+Qs-P00-Q00], [Qs-Q00 Qs-Q00]);
                end
            case SubEnh.None
        end

        F = get_F(Y, L, Ps, Qs, obj);
        nrm = norm(F);
        res = [res nrm];
        if k, kgmres = [kgmres gmresiter(end)]; end
        if ~silent, disp(['Iteration ' num2str(k) ': ' num2str(nrm) ' in ' num2str(gmresiter(end)) ' GMRES iterations (flag=' num2str(flag) ', relres=' num2str(relres) ')']), end
        if nrm < tol, break, end
        k = k + 1;

        if subenh.any && ~isempty(mp_c), mp_c.update_subenh(S, SP, SQ), end
        apply_jac_fun = get_apply_jac_fun(K, @subenh_prop_c, obj, N, DT);
        prec = get_prec(apply_jac_fun, K, @subenh_prop_c, mp_c, obj, N, DT, precinfo);

        [delta,flag,relres,gmresiter] = gmres(apply_jac_fun, -F, [], gmrestol, min(gmresmxit, numel(F)), prec);

        dY = reshape(delta(1:numel(delta)/2), d, []);
        dL = reshape(delta(numel(delta)/2+1:end), d, []);
        Y(:,2:1+size(dY,2)) = Y(:,2:1+size(dY,2)) + dY;
        L(:,2:1+size(dL,2)) = L(:,2:1+size(dL,2)) + dL;
    end

    function [P,Q] = subenh_prop_c(dy, dl, tstart, tend, obj, K, normalize)
        if subenh.any
            nn = tend / DT;
            
            comps = S' * [dy; dl];
            ylapprox = S * comps; yapprox = ylapprox(1:d); lapprox = ylapprox(d+1:end);
            [P,Q] = prop_c(dy - yapprox, dl - lapprox, tstart, tend, obj, K, normalize);
            P = P + SP * comps;
            Q = Q + SQ * comps;
            if ~normalize
                P = P - Pc00(:,nn) + P00(:,nn);
                Q = Q - Qc00(:,nn) + Q00(:,nn);
            end
        else
            [P,Q] = prop_c(dy, dl, tstart, tend, obj, K, normalize);
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

function prec = get_prec(A, K, prop_c, mp_c, obj, N, DT, precinfo)
    d = size(K, 1);

    if isempty(precinfo)
        prec = [];
        return
    end

    if precinfo.test
        prec = get_test_prec(A, d, N, obj, precinfo);
        return
    end

    switch precinfo.type
        case PrecType.Square
            prec = @(vec) square_prec(vec, K, prop_c, mp_c, obj, N, DT, precinfo);
        case PrecType.Triangular
            prec = @(vec) triangular_prec(vec, K, prop_c, mp_c, obj, N, DT, precinfo);
    end
end

function prec = get_test_prec(Afun, d, N, obj, precinfo)
    switch obj.type
        case ObjType.TerminalCost, n = d*N*2;
        case ObjType.Tracking, n = d*(N-1)*2;
    end

    A = sparse(n, n);
    for i=1:n
        v = zeros(n, 1);
        v(i) = 1;
        A(:,i) = Afun(v);
    end

    switch precinfo.type
        case PrecType.Square
            error 'Not implemented'
        case PrecType.Triangular
            prec = A;
            prec(n/2+1:end,1:n/2) = 0;
            prec(1:d,n/2-d+1:n/2) = precinfo.alpha * prec(d+1:2*d,1:d);
            prec(n-d+1:n,n/2+1:n/2+d) = precinfo.alpha' * prec(n/2+1:n/2+d,n/2+d+1:n/2+2*d);
    end
end

function res = square_prec(vec, K, prop_c, mp_c, obj, N, DT, precinfo)
    d = size(K, 1);
    switch obj.type
        case ObjType.Tracking
            M = N - 1;
        case ObjType.TerminalCost
            M = N;
    end

    % Step 1: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * precinfo.alpha^((m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * precinfo.alpha^((m-1)/M);
    end

    % Step 2: F
    vec = reshape(vec, d, 2*M);
    vec = [ifft(vec(:,1:M).'); ifft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)*sqrt(M);

    % Step 3: Solve systems
    D = zeros(M,1); D(2) = -precinfo.alpha^(1/M); D = M*ifft(D);
    for m=1:M
        if isempty(mp_c)
            [sol,~,~,~] = gmres(...
                @subsys, ...
                [vec((m-1)*d+1:m*d); vec(M*d+(m-1)*d+1:M*d+m*d)], ...
                [], 1e-10, 2*d...
            );
        else
            sol = [
                mp_c.I + D(m)*mp_c.Phi_f, mp_c.Psi_f;
                -mp_c.Psi_b, mp_c.I + D(m)'*mp_c.Phi_b;
            ] \ [mp_c.I*vec((m-1)*d+1:m*d); mp_c.I*vec(M*d+(m-1)*d+1:M*d+m*d)];
        end
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
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * (precinfo.alpha^(-(m-1)/M));
    end

    res = vec;

    function prd = subsys(v)
        prd = v;
        [P,~] = prop_c(D(m)*v(1:d), -v(d+1:end), m*DT, (m+1)*DT, obj, K, true);
        prd(1:d) = prd(1:d) + P;
        [~,Q] = prop_c(-v(1:d), D(m)'*v(d+1:end), m*DT, (m+1)*DT, obj, K, true);
        prd(d+1:end) = prd(d+1:end) + Q;
    end
end

function res = triangular_prec(vec, K, prop_c, mp_c, obj, N, DT, precinfo)
    d = size(K, 1);
    M = N;

    vec1 = vec(1:end/2);
    vec2 = vec(end/2+1:end);

    % PHASE 1: invert the bottom-right block
    %  Step 1: Gamma_alpha
    for m=1:M
        vec2((m-1)*d+1:m*d) = vec2((m-1)*d+1:m*d) * 1/(precinfo.alpha^((m-1)/M))';
    end
    %  Step 2: F
    vec2 = reshape(vec2, d, M);
    vec2 = ifft(vec2.').';
    vec2 = reshape(vec2, M*d, 1)*sqrt(M);
    %  Step 3: Solve systems
    D = zeros(M,1); D(2) = -precinfo.alpha^(1/M); D = M*ifft(D);
    vec2 = reshape(vec2, d, M);
    for m=1:M
        if isempty(mp_c)
            [sol,~,~,~] = gmres(...
                @subsys2, ...
                vec2(:,m), ...
                [], 1e-10, d...
            );
            vec2(:,m) = sol;
        else
            vec2(:,m) = (mp_c.I + D(m)'*mp_c.Phi_b)\(mp_c.I*vec2(:,m));
        end
    end
    vec2 = vec2(:);
    %  Step 4: F'
    vec2 = reshape(vec2, d, M);
    vec2 = fft(vec2.').';
    vec2 = reshape(vec2, M*d, 1)/sqrt(M);
    %  Step 5: Gamma_alpha
    for m=1:M
        vec2((m-1)*d+1:m*d) = vec2((m-1)*d+1:m*d) * (precinfo.alpha^((m-1)/M))';
    end
    
    % PHASE 2: invert the rest of the matrix
    %  Step 0: update with previous solution
    for m=1:M
        [P,~] = prop_c(zeros(d,1), vec2((m-1)*d+1:m*d), m*DT, (m+1)*DT, obj, K, true);
        vec1((m-1)*d+1:m*d) = vec1((m-1)*d+1:m*d) + P;
    end
    
    %  Step 1: Gamma_alpha
    for m=1:M
        vec1((m-1)*d+1:m*d) = vec1((m-1)*d+1:m*d) * precinfo.alpha^((m-1)/M);
    end
    %  Step 2: F
    vec1 = reshape(vec1, d, M);
    vec1 = ifft(vec1.').';
    vec1 = reshape(vec1, M*d, 1)*sqrt(M);
    %  Step 3: Solve systems
    D = zeros(M,1); D(2) = -precinfo.alpha^(1/M); D = M*ifft(D);
    vec1 = reshape(vec1, d, M);
    for m=1:M
        if isempty(mp_c)
            [sol,~,~,~] = gmres(...
                @subsys1, ...
                vec1(:,m), ...
                [], 1e-10, d...
            );
            vec1(:,m) = sol;
        else
            vec1(:,m) = (mp_c.I + D(m)*mp_c.Phi_f)\(mp_c.I*vec1(:,m));
        end
    end
    vec1 = vec1(:);
    %  Step 4: F'
    vec1 = reshape(vec1, d, M);
    vec1 = fft(vec1.').';
    vec1 = reshape(vec1, M*d, 1)/sqrt(M);
    %  Step 5: Gamma_alpha
    for m=1:M
        vec1((m-1)*d+1:m*d) = vec1((m-1)*d+1:m*d) * precinfo.alpha^(-(m-1)/M);
    end
    
    res = [vec1; vec2];

    function prd = subsys1(v)
        [P,~] = prop_c(D(m)*v, zeros(d,1), m*DT, (m+1)*DT, obj, K, true);
        prd = v + P;
    end

    function prd = subsys2(v)
        [~,Q] = prop_c(zeros(d,1), D(m)'*v, m*DT, (m+1)*DT, obj, K, true);
        prd = v + Q;
    end
end

function apply_jac_fun = get_apply_jac_fun(K, prop_c, obj, N, DT)
    switch obj.type
        case ObjType.Tracking, apply_jac_fun = @(delta) apply_jac_track(delta, K, prop_c, obj, N, DT);
        case ObjType.TerminalCost, apply_jac_fun = @(delta) apply_jac_tc(delta, K, prop_c, obj, N, DT);
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

%% ParaDiag
% Perform ParaDiag to calculate the optimal control u(t) of the equation
%   (1)  y'(t) = -Ky(t) + u(t),    y(0) = y0
% under a certain objective function.
%
% Parameters:
%  - K:         Matrix in the ODE (1) to control
%  - N:         There are `N+1` time discretisation points, at distance
%               `Tend/N` from each other
%  - Tend:      The end of the time interval [0, Tend] in which to solve
%               the problem
%  - y0:        The initial value in the ODE (1)
%  - obj:       Info about the objective function (see Obj class)
%  - precinfo:  Info about the preconditioner to use
%                 Default: No preconditioner
%  - showeigs:  Undocumented (used for creating figures)
%
% Return values:
%  - Y:      Discretised solution to (1)
%  - L:      Discretised adjoint
%  - kgmres: The number of GMRES iterations
%
function [Y,L,kgmres] = paradiag(K, N, Tend, y0, obj, precinfo, showeigs)
    if ~exist('precinfo', 'var') || isempty(precinfo), precinfo = []; end
    if ~exist('showeigs', 'var') || isempty(showeigs), showeigs = false; end

    dt = Tend / N;

    [A,b] = construct_system(K, N, dt, y0, obj);
    prec = get_prec(A, K, N, Tend, obj, precinfo);

    if showeigs
        assert(precinfo.test)
        [V,~]=eig(full(prec\A));cond(V)
        e = eig(full(prec\A)); scatter(real(e), imag(e), 60, 'o', 'filled')
    end

    if isempty(prec), maxiter = size(A, 1); else, maxiter = 25; end
    [sol,flag,relres,iter] = gmres(A, b, [], [], maxiter, prec);

    kgmres = iter(end);
    disp(['Solution of size-' num2str(size(A,1)) ' system found in ' num2str(kgmres) ' GMRES iterations (flag=' num2str(flag) ', relres=' num2str(relres) ')'])
    [Y,L] = postprocess(sol, y0, obj, K, dt);
end

function [Y,L] = postprocess(sol, y0, obj, K, dt)
    d = size(K, 1);

    switch obj.type
        case ObjType.Tracking
            Y = [y0 reshape(sol(1:numel(sol)/2), d, []) NaN(d,1)];
            Y(:,end) = (eye(d) + dt*K) \ Y(:,end-1);
            L = [NaN(d,1) reshape(sol(numel(sol)/2+1:end), d, []) zeros(d,1)];
            L(:,1) = (eye(d) + dt*K') \ (L(:,2) + dt/sqrt(obj.gamma)*(y0-obj.y_d(0)));
        case ObjType.TerminalCost
            Y = [y0 reshape(sol(1:numel(sol)/2), d, [])];
            L = [NaN(d,1) reshape(sol(numel(sol)/2+1:end), d, [])];
            L(:,1) = (eye(d) + dt*K') \ L(:,2);
    end
end

function prec = get_prec(A, K, N, Tend, obj, precinfo)
    d = size(K, 1);

    if isempty(precinfo)
        prec = [];
        return
    end

    if precinfo.test
        prec = get_test_prec(A, d, precinfo);
        return
    end

    switch precinfo.type
        case PrecType.Square
            prec = @(vec) square_prec(vec, K, N, d, Tend, obj, precinfo);
        case PrecType.Triangular
            prec = @(vec) triangular_prec(vec, K, N, d, Tend, obj, precinfo);
    end
end

function res = square_prec(vec, K, N, d, Tend, obj, precinfo)
    M = N - 1;
    dt = Tend / N;

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
    D = zeros(M,1); D(1:2) = [1 -1*precinfo.alpha^(1/M)]; D = M*ifft(D);
    if obj.type == ObjType.Tracking, gh = dt/sqrt(obj.gamma); else, gh = dt/obj.gamma; end
    v = [reshape(vec(1:M*d), d, []); reshape(vec(M*d+1:end), d, [])];
    for m=1:M
        sol = [
            D(m)*speye(d) + dt*K, gh*speye(d);
            -gh*speye(d)*(obj.type == ObjType.Tracking), D(m)'*speye(d) + dt*K';
        ] \ v(:,m);
        v(:,m) = sol;
    end
    v1 = v(1:d,:); v2 = v(d+1:end,:); vec = [v1(:); v2(:)];
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
end

function res = triangular_prec(vec, K, N, d, Tend, obj, precinfo)
    M = N;
    dt = Tend / N;
    
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
    D = zeros(M,1); D(1:2) = [1 -1*precinfo.alpha^(1/M)]; D = M*ifft(D);
    vec2 = reshape(vec2, d, M);
    for m=1:M
        vec2(:,m) = (D(m)'*speye(d)+dt*K')\vec2(:,m);
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
    vec1 = vec1 - dt/obj.gamma*vec2;
    %  Step 1: Gamma_alpha
    for m=1:M
        vec1((m-1)*d+1:m*d) = vec1((m-1)*d+1:m*d) * precinfo.alpha^((m-1)/M);
    end
    %  Step 2: F
    vec1 = reshape(vec1, d, M);
    vec1 = ifft(vec1.').';
    vec1 = reshape(vec1, M*d, 1)*sqrt(M);
    %  Step 3: Solve systems
    D = zeros(M,1); D(1:2) = [1 -1*precinfo.alpha^(1/M)]; D = M*ifft(D);
    vec1 = reshape(vec1, d, M);
    for m=1:M
        vec1(:,m) = (D(m)*speye(d)+dt*K)\vec1(:,m);
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
end

function prec = get_test_prec(A, d, precinfo)
    n = size(A, 1);
    prec = A;

    prec(1:d,n/2-d+1:n/2) = precinfo.alpha * prec(d+1:2*d,1:d);
    prec(n-d+1:n,n/2+1:n/2+d) = precinfo.alpha' * prec(n/2+1:n/2+d,n/2+d+1:n/2+2*d);

    switch precinfo.type
        case PrecType.Triangular
            prec(end-d+1:end,n/2-d+1:n/2) = 0;
        case PrecType.Square
    end
end

function [A,b] = construct_system(K, N, dt, y0, obj)
    d = size(K, 1);

    switch obj.type
        case ObjType.Tracking, M = N-1;
        case ObjType.TerminalCost, M = N;
    end

    K = sparse(K);
    A = blkdiag(kron(speye(M), speye(d) + dt*K), kron(speye(M), speye(d) + dt*K'));
    A = A - [
        spdiags(ones((M-1)*d, 1), -d, M*d, M*d), sparse(M*d, M*d);
        sparse(M*d, M*d), spdiags(ones((M-1)*d, 1), -d, M*d, M*d)';
    ];
    b = sparse(2*M*d,1);
    b(1:d) = y0;

    switch obj.type
        case ObjType.Tracking
            A = A + [
                sparse(M*d, M*d), dt/sqrt(obj.gamma)*speye(M*d);
                -dt/sqrt(obj.gamma)*speye(M*d), sparse(M*d, M*d);
            ];
            for m=1:M
                b(M*d+(m-1)*d+1:M*d+m*d) = -dt/sqrt(obj.gamma)*obj.y_d(m*dt);
            end
        case ObjType.TerminalCost
            A = A + [
                sparse(M*d, M*d), dt/obj.gamma*speye(M*d);
                sparse(M*d, 2*M*d);
            ];
            A(end-d+1:end,M*d-d+1:M*d) = -speye(d) - dt*K';
            b(end-d+1:end) = (-speye(d) - dt*K')*obj.y_T;
    end
end

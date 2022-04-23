%% ParaDiag
%
function [Y,L] = paradiag(K, N, Tend, y0, obj, prectype)
    if ~exist('prectype','var') || isempty(prectype), prectype = []; end
    if obj.type ~= ObjType.TerminalCost, assert(abs(abs(prectype.alpha)-1) < 1e-12), end

    dt = Tend / N;

    [A,b] = construct_system(K, N, dt, y0, obj);
    prec = get_prec(A, K, N, Tend, obj, prectype);
    
%     n = size(A,1);
%     prc = zeros(n,n);
%     for i=1:n
% %         disp([num2str(i) ' / ' num2str(n)])
%         v = zeros(n,1); v(i) = 1;
%         prc(:,i) = prec(v);
%     end
%     prectype.test = true;
%     prec = get_prec(A, K, N, Tend, obj, prectype);
% %     imshow(abs(inv(prc))), figure, imshow(abs(full(prec)))
% %     imshow(abs(prc*prec)),pause
%     imshow(abs(inv(prc)-prec)), norm(abs(inv(prc)-prec)), figure, imshow(abs(prc*prec)), pause
% %     condest(prc\A)

    [sol, flag, relres, iter] = gmres(A, b, [], [], size(A,1), prec);

    disp(['Solution of size-' num2str(size(A,1)) ' system found in ' num2str(iter(end)) ' GMRES iterations (flag=' num2str(flag) ', relres=' num2str(relres) ')'])
    [Y,L] = postprocess(sol, y0, obj, K, dt);
end

function [Y,L] = postprocess(sol, y0, obj, K, dt)
    d = size(K, 1);

    switch obj.type
        case ObjType.Tracking
            Y = [y0 reshape(sol(1:numel(sol)/2), d, []) NaN(d,1)];
            Y(:,end) = (eye(d) + dt*K) \ Y(:,end-1);
            L = [NaN(d,1) reshape(sol(numel(sol)/2+1:end), d, []) zeros(d,1)];
            L(:,1) = (eye(d) + dt*K) \ L(:,2);
        case ObjType.TerminalCost
            Y = [y0 reshape(sol(1:numel(sol)/2), d, [])];
            L = [NaN(d,1) reshape(sol(numel(sol)/2+1:end), d, [])];
            L(:,1) = (eye(d) + dt*K) \ L(:,2);
        case ObjType.TerminalCostMod
            Y = [y0 reshape(sol(1:numel(sol)/2), d, []) NaN(d,1)];
            L = [NaN(d,1) reshape(sol(numel(sol)/2+1:end), d, []) -obj.y_T];
            L(:,1) = (eye(d) + dt*K + dt/obj.gamma*eye(d)) \ (L(:,2) + (- 2*dt*K - dt/obj.gamma*eye(d))*y0);
            Y(:,end) = (eye(d) + dt*K + dt/obj.gamma*eye(d)) \ (Y(:,end-1) - dt/obj.gamma*L(:,end));
    end
end

function prec = get_prec(A, K, N, Tend, obj, prectype)
    d = size(K, 1);

    if isempty(prectype)
        prec = [];
        return
    end

    if prectype.test
        prec = get_test_prec(A, d, obj, prectype);
        return
    end

    switch obj.type
        case ObjType.Tracking
            prec = @(vec) tracking_prec(vec, K, N, d, Tend, obj, prectype);
        case ObjType.TerminalCost
            prec = @(vec) tc_prec(vec, K, N, d, Tend, obj, prectype);
        case ObjType.TerminalCostMod
            prec = @(vec) tc_mod_prec(vec, K, N, d, Tend, obj, prectype);
    end
end

function res = tracking_prec(vec, K, N, d, Tend, obj, prectype)
    M = N - 1;
    dt = Tend / N;

    % Step 1: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * prectype.alpha^((m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * 1/(prectype.alpha^((m-1)/M))';
    end
    % Step 2: F
    vec = reshape(vec, d, 2*M);
    vec = [ifft(vec(:,1:M).'); ifft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)*sqrt(M);
    % Step 3: Solve systems
    D = zeros(M,1); D(1:2) = [1 -1*prectype.alpha^(1/M)]; D = M*ifft(D);
    for m=1:M
        sol = [
            D(m)*eye(d) + dt*K, dt/sqrt(obj.gamma)*eye(d)*(1-prectype.diag);
            -dt/sqrt(obj.gamma)*eye(d)*(1-prectype.diag), D(m)'*eye(d) + dt*K;
        ] \ [vec((m-1)*d+1:m*d); vec(M*d+(m-1)*d+1:M*d+m*d)];
        vec((m-1)*d+1:m*d) = sol(1:d);
        vec(M*d+(m-1)*d+1:M*d+m*d) = sol(d+1:end);
    end
    % Step 4: F'
    vec = reshape(vec, d, 2*M);
    vec = [fft(vec(:,1:M).'); fft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)/sqrt(M);
    % Step 5: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * prectype.alpha^(-(m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * (prectype.alpha^((m-1)/M))';
    end

    res = vec;
end

function res = tc_prec(vec, K, N, d, Tend, obj, prectype)
    M = N;
    dt = Tend / N;

    % Step 1: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * prectype.alpha^((m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * 1/(prectype.alpha^((m-1)/M))';
    end
    % Step 2: F
    vec = reshape(vec, d, 2*M);
    vec = [ifft(vec(:,1:M).'); ifft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)*sqrt(M);
    % Step 3: Solve systems
    D = zeros(M,1); D(1:2) = [1 -1*prectype.alpha^(1/M)]; D = M*ifft(D);
    for m=1:M
        sol = [
            D(m)*eye(d) + dt*K, dt/obj.gamma*eye(d)*(1-prectype.diag);
            zeros(d), D(m)'*eye(d) + dt*K;
        ] \ [vec((m-1)*d+1:m*d); vec(M*d+(m-1)*d+1:M*d+m*d)];
        vec((m-1)*d+1:m*d) = sol(1:d);
        vec(M*d+(m-1)*d+1:M*d+m*d) = sol(d+1:end);
    end
    % Step 4: F'
    vec = reshape(vec, d, 2*M);
    vec = [fft(vec(:,1:M).'); fft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)/sqrt(M);
    % Step 5: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * prectype.alpha^(-(m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * (prectype.alpha^((m-1)/M))';
    end

    res = vec;
end

function res = tc_mod_prec(vec, K, N, d, Tend, obj, prectype)
    M = N - 1;
    dt = Tend / N;

    % Step 1: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * prectype.alpha^((m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * 1/(prectype.alpha^((m-1)/M))';
    end
    % Step 2: F
    vec = reshape(vec, d, 2*M);
    vec = [ifft(vec(:,1:M).'); ifft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)*sqrt(M);
    % Step 3: Solve systems
    D = zeros(M,1); D(1:2) = [1 -1*prectype.alpha^(1/M)]; D = M*ifft(D);
    for m=1:M
        sol = [
            D(m)*eye(d) + dt*K + dt/obj.gamma*eye(d), dt/obj.gamma*eye(d)*(1-prectype.diag);
            dt/obj.gamma*speye(d)+2*dt*K, D(m)'*eye(d) + dt*K + dt/obj.gamma*eye(d);
        ] \ [vec((m-1)*d+1:m*d); vec(M*d+(m-1)*d+1:M*d+m*d)];
        vec((m-1)*d+1:m*d) = sol(1:d);
        vec(M*d+(m-1)*d+1:M*d+m*d) = sol(d+1:end);
    end
    % Step 4: F'
    vec = reshape(vec, d, 2*M);
    vec = [fft(vec(:,1:M).'); fft(vec(:,M+1:end).')].';
    vec = reshape(vec, 2*M*d, 1)/sqrt(M);
    % Step 5: Gamma_alpha
    for m=1:M
        vec((m-1)*d+1:m*d) = vec((m-1)*d+1:m*d) * prectype.alpha^(-(m-1)/M);
        vec(M*d+(m-1)*d+1:M*d+m*d) = vec(M*d+(m-1)*d+1:M*d+m*d) * (prectype.alpha^((m-1)/M))';
    end

    res = vec;
end

function prec = get_test_prec(A, d, obj, prectype)
    n = size(A, 1);

    prec = sparse(n, n);
    for i=1:n
        v = zeros(n, 1);
        v(i) = 1;
        prec(:,i) = A*v;
    end

    switch obj.type
        case {ObjType.Tracking, ObjType.TerminalCostMod}
            prec(1:d,n/2-d+1:n/2) = prectype.alpha * prec(d+1:2*d,1:d);
            prec(n-d+1:n,n/2+1:n/2+d) = prectype.alpha' * prec(n/2+1:n/2+d,n/2+d+1:n/2+2*d);
            if prectype.diag
                prec(n/2+1:end,1:n/2) = 0; prec(1:n/2,n/2+1:end) = 0;
            end
        case ObjType.TerminalCost
            prec = A;
            prec(end-d+1:end,n/2-d+1:n/2) = 0;
            prec(1:d,n/2-d+1:n/2) = prectype.alpha * prec(d+1:2*d,1:d);
            prec(end-d+1:end,n/2+1:n/2+d) = prectype.alpha' * prec(n/2+1:n/2+d,n/2+d+1:n/2+2*d);
            if prectype.diag
                prec(n/2+1:end,1:n/2) = 0; prec(1:n/2,n/2+1:end) = 0;
            end

%             mat111 = [speye(n/2), sparse(n/2, n/2); speye(n/2), speye(n/2)];
%             mat1m11 = [speye(n/2), sparse(n/2, n/2); -speye(n/2), speye(n/2)];
%             prec = A;
%             prec(1:d,n/2-2*d+1:n/2-d) = prectype.alpha*prec(d+1:2*d,1:d);
%             prec(end-2*d+1:end-d,n/2+1:n/2+d) = prectype.alpha*prec(n/2+1:n/2+d,n/2+d+1:n/2+2*d);
%             prec = mat111*prec*mat111;
% 
% %             prec(n/2+1:end,1:n/2) = prec(n/2+1:end,1:n/2) .* kron(eye(n/2/d), ones(d));
% %             prec(n/2+1:end,1:n/2) = prec(n/2+1:end,1:n/2) - 2*speye(n/2);
% 
%             prec(end-d+1:end,n/2-d+1:n/2) = 0;
%             prec(end-d+1:end,end-d+1:end) = eye(d);
% 
%             prec = mat1m11*prec*mat1m11;
    end
end

function [A,b] = construct_system(K, N, dt, y0, obj)
    d = size(K, 1);

    switch obj.type
        case {ObjType.Tracking, ObjType.TerminalCostMod}
            M = N-1;
        case ObjType.TerminalCost
            M = N;
    end

    K = sparse(K);
    A = kron(speye(2*M), speye(d) + dt*K);
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
            A(end-d+1:end,M*d-d+1:M*d) = -speye(d) - dt*K;
            b(end-d+1:end) = (-speye(d) - dt*K)*obj.y_T;
        case ObjType.TerminalCostMod
            A = A + [
                kron(speye(M), dt/obj.gamma*speye(d)), dt/obj.gamma*speye(M*d);
                kron(speye(M), dt/obj.gamma*speye(d)+2*dt*K), kron(speye(M), dt/obj.gamma*speye(d))
            ];
            b(end-d+1:end) = -obj.y_T;
    end
end

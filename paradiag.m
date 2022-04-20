%% ParaDiag
%
function [Y,L] = paradiag(K, N, Tend, y0, obj)
    dt = Tend / N;

    [A,b] = construct_system(K, N, dt, y0, obj);
    sol = full(A\b);
    [Y,L] = postprocess(sol, y0, obj, K, dt);
end

function [Y,L] = postprocess(sol, y0, obj, K, dt)
    d = size(K, 1);

    switch obj.type
        case ObjType.Tracking
            Y = [y0 reshape(sol(1:numel(sol)/2), d, []) NaN(d,1)];
            Y(:,end) = (eye(d) + dt*K)\Y(:,end-1);
            L = [NaN(d,1) reshape(sol(numel(sol)/2+1:end), d, []) zeros(d,1)];
            L(:,1) = (eye(d) + dt*K)\L(:,2);
        case ObjType.TerminalCost
            Y = [y0 reshape(sol(1:numel(sol)/2), d, [])];
            L = [NaN(d,1) reshape(sol(numel(sol)/2+1:end), d, [])];
            L(:,1) = (eye(d) + dt*K)\L(:,2);
    end
end

function [A,b] = construct_system(K, N, dt, y0, obj)
    d = size(K, 1);

    switch obj.type
        case ObjType.Tracking
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
    end
end

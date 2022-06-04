%% Tracking implicit Euler propagator
% This uses implicit Euler to solve the ParaOpt sub-problems for tracking
% problems.
%
function [yend,l0] = prop_ie_track(steps, y0, lend, Tstart, Tend, obj, K, normalize)
    d = size(K,1);
    fll = 2*(steps+1)*d;
    half = fll / 2;

    dt = (Tend - Tstart) / steps;
    IpKdt = speye(d) + K*dt;
    IpKtdt = speye(d) + K'*dt;

    b = sparse(fll,1);
    b(1:d) = y0;
    b(end-d+1:end) = lend;

    if ~normalize
        for i=1:steps
            idx = (i-1)*d+1:i*d;
            b(idx+half) = eye(d)*dt/sqrt(obj.gamma)*obj.y_d(Tstart + i*dt);
        end
    end
    
    M = blkdiag(kron(speye(half/d), -IpKdt), kron(speye(half/d), -IpKtdt)) + [
        sparse(d,fll);
        [
            speye(steps*d), sparse(steps*d,2*d), -dt/sqrt(obj.gamma)*speye(steps*d);
            dt/sqrt(obj.gamma)*speye(steps*d), sparse(steps*d,2*d), speye(steps*d);
        ];
        sparse(d,fll);
    ];
    M(1:d,1:d) = speye(d);
    M(end-d+1:end,end-d+1:end) = speye(d);

    solved = M\b;
    yend = solved(half-d+1:half);
    l0 = solved(half+1:half+d);
end

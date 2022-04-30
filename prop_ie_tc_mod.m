% TODO: Construct the sparse matrices faster
function [yend,l0] = prop_ie_tc_mod(steps, y0, lend, Tstart, Tend, obj, K, ~)
    d = size(K,1);
    fll = 2*(steps+1)*d;
    half = fll / 2;

    dt = (Tend - Tstart) / steps;
    L1 = speye(d) + K*dt + speye(d)*dt/obj.gamma;
    L2 = -2*K*dt - speye(d)*dt/obj.gamma;

    M = sparse(fll,fll);
    b = sparse(fll,1);
    M(1:d,1:d) = speye(d); b(1:d) = y0;
    M(end-d+1:end,end-d+1:end) = speye(d); b(end-d+1:end) = lend;

    for i=1:steps
        idx = (i-1)*d+1:i*d;
        M(idx+d,idx) = eye(d);
        M(idx+d,idx+d) = -L1;
        M(idx+d,idx+half+d) = -dt/obj.gamma*eye(d);
    end

    for i=1:steps
        idx = (i-1)*d+1:i*d;
        M(idx+half,idx) = L2;
        M(idx+half,idx+half) = -L1;
        M(idx+half,idx+half+d) = eye(d);
    end

    solved = M\b;
    yend = solved(half-d+1:half);
    l0 = solved(half+1:half+d);
end

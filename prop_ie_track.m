% TODO: Construct the sparse matrices faster
function [yend,l0] = prop_ie_track(steps, y0, lend, Tstart, Tend, obj, K, deriv)
    d = size(K,1);
    fll = 2*(steps+1)*d;
    half = fll / 2;

    dt = (Tend - Tstart) / steps;
    IpKdt = speye(d) + K*dt;

    M = sparse(fll,fll);
    b = sparse(fll,1);
    M(1:d,1:d) = speye(d); b(1:d) = y0;
    M(end-d+1:end,end-d+1:end) = speye(d); b(end-d+1:end) = lend;

    for i=1:steps
        idx = (i-1)*d+1:i*d;
        M(idx+d,idx) = eye(d);
        M(idx+d,idx+d) = -IpKdt;
        M(idx+d,idx+half+d) = -dt/sqrt(obj.gamma)*eye(d);
    end

    for i=1:steps
        idx = (i-1)*d+1:i*d;
        M(idx+half,idx) = eye(d)*dt/sqrt(obj.gamma);
        M(idx+half,idx+half) = -IpKdt;
        M(idx+half,idx+half+d) = eye(d);
        if ~deriv
            b(idx+half) = eye(d)*dt/sqrt(obj.gamma)*obj.y_d(Tstart + i*dt);
        end
    end

    solved = M\b;
    yend = solved(half-d+1:half);
    l0 = solved(half+1:half+d);
end

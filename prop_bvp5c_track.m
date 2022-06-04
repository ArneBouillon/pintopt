% TODO: y_d
function [yend,l0] = prop_bvp5c_track(y0, lend, Tstart, Tend, obj, K, normalize, tol, pts)
    if ~exist('tol', 'var') || isempty(tol), tol = 1e-6; end
    if ~exist('pts', 'var') || isempty(pts), pts = 5; end

    d = size(K,1);

    f = @(t,yl) [
        -K*yl(1:d) - yl(d+1:end)/sqrt(obj.gamma);
        (obj.y_d(t)*(1-normalize)-yl(1:d))/sqrt(obj.gamma) + K'*yl(d+1:end);
    ];
    fjac = sparse([-K -eye(d)/sqrt(obj.gamma); -eye(d)/sqrt(obj.gamma) K]);
    bounds = @(yla, ylb) [
        yla(1:d) - y0;
        ylb(d+1:end) - lend;
    ];
    bcjac = {[speye(d) sparse(d,d); sparse(d,2*d)], [sparse(d,2*d); sparse(d,d) speye(d)]};

    solinit = struct;
    solinit.x = linspace(Tstart, Tend, pts);
    solinit.y = [y0*ones(1, pts); lend*ones(1, pts)];
    sol = bvp5c(f, bounds, solinit, bvpset("RelTol", tol, "FJacobian", fjac, "BCJacobian", bcjac));
    yend = sol.y(1:d,end);
    l0 = sol.y(d+1:2*d,1);
end

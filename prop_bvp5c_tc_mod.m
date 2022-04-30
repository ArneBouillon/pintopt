function [yend,l0] = prop_bvp5c_tc_mod(tol, pts, y0, lend, Tstart, Tend, obj, K, ~)
    d = size(K,1);

    f = @(t,yl) [
        (-K-eye(d)/obj.gamma)*yl(1:d) - yl(d+1:end)/obj.gamma;
        (2*K+eye(d)/obj.gamma)*yl(1:d) + (K+eye(d)/obj.gamma)*yl(d+1:end);
    ];
    fjac = sparse([-K-eye(d)/obj.gamma, -eye(d)/obj.gamma; 2*K+eye(d)/obj.gamma, K+eye(d)/obj.gamma]);
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

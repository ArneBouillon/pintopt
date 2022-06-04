%% Terminal-cost BVP5C propagator
% This uses MATLAB's built-in BVP5C boundary-value--problem solver to solve
% the ParaOpt sub-problems for terminal-cost problems.
%
function [yend,l0] = prop_bvp5c_tc(y0, lend, Tstart, Tend, obj, K, ~, tol, pts)
    if ~exist('tol', 'var') || isempty(tol), tol = 1e-6; end
    if ~exist('pts', 'var') || isempty(pts), pts = 5; end

    d = size(K,1);

    f = @(t,yl) [
        -K*yl(1:d) - yl(d+1:end)/obj.gamma;
        K'*yl(d+1:end);
    ];
    fjac = sparse([-K -eye(d)/obj.gamma; zeros(d) K]);
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

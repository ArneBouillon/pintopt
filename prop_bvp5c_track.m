function [yend,l0] = prop_bvp5c_track(tol, pts, Tstart, Tend, obj, A, y0, lend, deriv)
    d = size(A,1);

    f = @(t,yl) [
        A*yl(1:d) - yl(d+1:end)/sqrt(obj.gamma);
        -yl(1:d)/sqrt(obj.gamma) - A*yl(d+1:end) + (obj.y_d(t)/sqrt(obj.gamma))*(~deriv);
    ];
    bounds = @(yla, ylb) [
        yla(1:d) - y0;
        ylb(d+1:end) - lend;
    ];

    solinit = struct;
    solinit.x = linspace(Tstart, Tend, pts);
    solinit.y = [y0*ones(1, pts); lend*ones(1, pts)];
    sol = bvp5c(f, bounds, solinit, bvpset("RelTol", tol));
    yend = sol.y(1:d,end);
    l0 = sol.y(d+1:2*d,1);
end

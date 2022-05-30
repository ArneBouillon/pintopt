function [yend,l0] = prop_trap1_track(y0, lend, Tstart, Tend, obj, K, normalize)
    d = size(K,1);
    dt = Tend - Tstart;
    Ix = speye(d);
    Ox = sparse(d, d);

    M = [
        Ix Ox Ox Ox;
        Ix-dt*K/2, -(Ix+dt*K/2), -dt/2/sqrt(obj.gamma)*Ix, -dt/2/sqrt(obj.gamma)*Ix;
        dt/2/sqrt(obj.gamma)*Ix, dt/2/sqrt(obj.gamma)*Ix, -(Ix+dt*K'/2), Ix-dt*K'/2
        Ox Ox Ox Ix;
    ];
    b = [
        y0;
        zeros(d,1);
        (1-normalize) * (dt/2/sqrt(obj.gamma)*(obj.y_d(Tstart) + obj.y_d(Tend)));
        lend;
    ];

    solved = M\b;
    yend = solved(d+1:2*d);
    l0 = solved(2*d+1:3*d);
end

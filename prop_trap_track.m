function [yend,l0] = prop_trap_track(steps, y0, lend, Tstart, Tend, obj, K, normalize)
    d = size(K,1);
    fll = 2*(steps+1)*d;
    half = fll / 2;

    dt = (Tend - Tstart) / steps;
    Ix = speye(d);
    Ox = sparse(d, d);
    
    M = blkdiag(kron(speye(half/d), -(Ix+dt*K/2)), kron(speye(half/d), -(Ix+dt*K'/2)));
    M(1:d,1:d) = speye(d); M(end-d+1:end,end-d+1:end) = speye(d);
    B = kron(diag(sparse(ones(half/d-1,1)),-1), Ix-dt*K/2);
    M = M + blkdiag(B, B');
    M = M + [
        sparse(d,fll);
        [
            sparse(steps*d,steps*d+2*d), -dt/sqrt(obj.gamma)/2*speye(steps*d);
            dt/sqrt(obj.gamma)/2*speye(steps*d), sparse(steps*d,steps*d+2*d);
        ];
        sparse(d,fll);
    ] + [
        sparse(d,fll);
        [
            sparse(steps*d,steps*d+d), -dt/sqrt(obj.gamma)/2*speye(steps*d), sparse(steps*d,d);
            sparse(steps*d,d), dt/sqrt(obj.gamma)/2*speye(steps*d), sparse(steps*d,steps*d+d)
        ];
        sparse(d,fll);
    ];

    MM = [
        Ix Ox Ox Ox;
        Ix-dt*K/2, -(Ix+dt*K/2), -dt/2/sqrt(obj.gamma)*Ix, -dt/2/sqrt(obj.gamma)*Ix;
        dt/2/sqrt(obj.gamma)*Ix, dt/2/sqrt(obj.gamma)*Ix, -(Ix+dt*K'/2), Ix-dt*K'/2
        Ox Ox Ox Ix;
    ];

%     spy(M), figure, spy(MM), figure, spy(M-MM), pause
%     b = [
%         y0;
%         zeros(d,1);
%         (1-normalize) * (dt/2/sqrt(obj.gamma)*(obj.y_d(Tstart) + obj.y_d(Tend)));
%         lend;
%     ];
    b = [
        y0;
        zeros(steps*d,1);
        zeros(steps*d,1);
        lend;
    ];
    if ~normalize
        for i=1:steps
            b((steps+i)*d+1:(steps+i+1)*d) = dt/2/sqrt(obj.gamma)*(obj.y_d(Tstart + (i-1)*dt) + obj.y_d(Tstart + i*dt));
        end
    end

    solved = M\b;
    yend = solved(d+1:2*d);
    l0 = solved(2*d+1:3*d);
end

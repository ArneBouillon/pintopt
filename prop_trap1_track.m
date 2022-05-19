function [yend,l0] = prop_trap1_track(y0, lend, Tstart, Tend, obj, K, deriv)
    d = size(K,1);
    dt = Tend - Tstart;
    Ix = speye(d);
    Ox = sparse(d, d);
    
    if deriv
        sh = K*dt; gh = Ix*dt/sqrt(obj.gamma);
        inv = (Ix+sh/2+(gh/2)^2/(Ix+sh/2));
        phi = inv\(Ix-sh/2-(gh/2)^2/(Ix+sh/2));
        psi = inv\(gh/2*(Ix+Ix/(Ix+sh/2)*(Ix-sh/2)));
        yend = phi*y0-psi*lend;
        l0 = psi*y0+phi*lend;
        return
    end

    M = [
        Ix Ox Ox Ox;
        Ix-dt*K/2, -(Ix+dt*K/2), -dt/2/sqrt(obj.gamma)*Ix, -dt/2/sqrt(obj.gamma)*Ix;
        dt/2/sqrt(obj.gamma)*Ix, dt/2/sqrt(obj.gamma)*Ix, -(Ix+dt*K/2), Ix-dt*K/2
        Ox Ox Ox Ix;
    ];
    b = [
        y0;
        zeros(d,1);
        (1-deriv) * (dt/2/sqrt(obj.gamma)*(obj.y_d(Tstart) + obj.y_d(Tend)));
        lend;
    ];

%     IpKdt = speye(d) + K*dt;
% 
%     b = sparse(fll,1);
%     b(1:d) = y0;
%     b(end-d+1:end) = lend;
% 
%     if ~deriv
%         for i=1:steps
%             idx = (i-1)*d+1:i*d;
%             b(idx+half) = eye(d)*dt/sqrt(obj.gamma)*obj.y_d(Tstart + i*dt);
%         end
%     end
%     
%     M = kron(speye(fll/d), -IpKdt) + [
%         sparse(d,fll);
%         [
%             speye(steps*d), sparse(steps*d,2*d), -dt/sqrt(obj.gamma)*speye(steps*d);
%             dt/sqrt(obj.gamma)*speye(steps*d), sparse(steps*d,2*d), speye(steps*d);
%         ];
%         sparse(d,fll);
%     ];
%     M(1:d,1:d) = speye(d);
%     M(end-d+1:end,end-d+1:end) = speye(d);

    solved = M\b;
    yend = solved(d+1:2*d);
    l0 = solved(2*d+1:3*d);
end

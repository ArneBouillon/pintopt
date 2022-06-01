% TODO: Transform this to [Gander]'s discretisation instead of our own?
function [yend,l0] = prop_ie_tc(steps, y0, lend, Tstart, Tend, obj, K, ~, fdto)
    if ~exist('fdto','var') || isempty(fdto)
        fdto = false;
    end

    d = size(K,1);
    dt = (Tend - Tstart) / steps;
    IpKdt = speye(d) + K*dt;
    IpKtdt = speye(d) + K'*dt;
    
    ls = NaN(d, steps+1);
    ls(:,end) = lend;
    for i=steps:-1:1
        ls(:,i) = IpKtdt\ls(:,i+1);
    end
    l0 = ls(:,1);
    
    ys = NaN(d, steps+1);
    ys(:,1) = y0;
    
    for i=2:steps+1
        if fdto
            ys(:,i) = IpKdt\(ys(:,i-1) - dt/obj.gamma*(IpKtdt\ls(:,i)));
        else
            ys(:,i) = IpKdt\(ys(:,i-1) - dt/obj.gamma*ls(:,i));
        end
    end
    yend = ys(:,end);
end

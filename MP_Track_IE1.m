classdef MP_Track_IE1 < MultProp
    properties
        K
        d
        obj
        dt
    end

    methods
        function mp = MP_Track_IE1(K, obj, dt)
            mp.K = K;
            mp.d = size(K, 1);
            mp.obj = obj;
            mp.dt = dt;
        end

        function I = I(self)
            I = eye(self.d) + self.dt*self.K;
        end

        function Phi = Phi(self)
            Phi = eye(self.d);
        end

        function Psi = Psi(self)
            Psi = self.dt / sqrt(self.obj.gamma) * eye(self.d);
        end
    end
end

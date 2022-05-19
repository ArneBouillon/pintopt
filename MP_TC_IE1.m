classdef MP_TC_IE1 < handle
    properties
        d
        I_
        Phi_
        Psi_
        I
        Phi_f
        Phi_b
        Psi_f
        Psi_b
    end

    methods
        function mp = MP_TC_IE1(K, obj, dt)
            d = size(K, 1);
            mp.d = d;

            mp.I_ = speye(d) + dt*K;
            mp.Phi_ = speye(d);
            mp.Psi_ = dt/obj.gamma*speye(d);

            mp.I = mp.I_;
            mp.Phi_f = mp.Phi_;
            mp.Phi_b = mp.Phi_;
            mp.Psi_f = mp.Psi_;
            mp.Psi_b = sparse(d,d);
        end

        function update_subenh(self, S, SP, SQ)
            nonproj = eye(2*self.d) - S*S';
            Pproj = SP*S';
            Qproj = SQ*S';
            self.Phi_f = self.Phi_ * nonproj(1:self.d,1:self.d) ...
                - self.Psi_ * nonproj(self.d+1:end,1:self.d) ...
                + self.I * Pproj(:,1:self.d);
            self.Psi_f = self.Psi_ * nonproj(self.d+1:end,self.d+1:end) ...
                - self.Phi_ * nonproj(1:self.d,self.d+1:end) ...
                - self.I * Pproj(:,self.d+1:end);
            self.Psi_b = self.Phi_ * nonproj(self.d+1:end,1:self.d) ...
                + self.I * Qproj(:,1:self.d);
            self.Phi_b = self.Phi_ * nonproj(self.d+1:end,self.d+1:end) ...
                + self.I * Qproj(:,self.d+1:end);
        end
    end
end

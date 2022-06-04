%% Multiplied propagator for single-step tracking implicit Euler
% See Appendix C of the thesis for information about MPs
%
classdef MP_Track_IE1 < handle
    properties
        d

        Xi_f_
        Xi_b_
        XiPhi_f_
        XiPhi_b_
        XiPsi_f_
        XiPsi_b_

        Xi_f
        Xi_b
        XiPhi_f
        XiPhi_b
        XiPsi_f
        XiPsi_b
    end

    methods
        function mp = MP_Track_IE1(K, obj, dt)
            d = size(K, 1);
            mp.d = d;

            mp.Xi_f_ = speye(d) + dt*K;
            mp.Xi_b_ = speye(d) + dt*K';
            mp.XiPhi_f_ = speye(d);
            mp.XiPhi_b_ = speye(d);
            mp.XiPsi_f_ = dt/sqrt(obj.gamma)*speye(d);
            mp.XiPsi_b_ = dt/sqrt(obj.gamma)*speye(d);
            
            mp.Xi_f = mp.Xi_f_;
            mp.Xi_b = mp.Xi_b_;
            mp.XiPhi_f = mp.XiPhi_f_;
            mp.XiPhi_b = mp.XiPhi_b_;
            mp.XiPsi_f = mp.XiPsi_f_;
            mp.XiPsi_b = mp.XiPsi_b_;
        end

        function update_subenh(self, S, SP, SQ)
            nonproj = eye(2*self.d) - S*S';
            Pproj = SP*S';
            Qproj = SQ*S';
            
            self.XiPhi_f = self.XiPhi_f_ * nonproj(1:self.d,1:self.d) ...
                - self.XiPsi_f_ * nonproj(self.d+1:end,1:self.d) ...
                + self.Xi_f_ * Pproj(:,1:self.d);
            self.XiPsi_f = self.XiPsi_f_ * nonproj(self.d+1:end,self.d+1:end) ...
                - self.XiPhi_f_ * nonproj(1:self.d,self.d+1:end) ...
                - self.Xi_f_ * Pproj(:,self.d+1:end);
            self.XiPsi_b = self.XiPsi_b_ * nonproj(1:self.d,1:self.d) ...
                + self.XiPhi_b_ * nonproj(self.d+1:end,1:self.d) ...
                + self.Xi_b_ * Qproj(:,1:self.d);
            self.XiPhi_b = self.XiPhi_b_ * nonproj(self.d+1:end,self.d+1:end) ...
                + self.XiPsi_b_ * nonproj(1:self.d,self.d+1:end) ...
                + self.Xi_b_ * Qproj(:,self.d+1:end);
        end
    end
end

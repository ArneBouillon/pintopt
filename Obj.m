%% Objective function class
% There are different types of objective function.
%  - Terminal cost: (1) J_T(y,u) = 1/2*||y(T)-y_T||^2 + gamma/2 * int_0^T(||u||^2)
%  - Tracking:      (2) J_d(y,u) = 1/2*int_0^T(||y-y_d||^2 + gamma*||u||^2)
%
% Parameters:
%  - type:  The type of objective function
%  - gamma: The parameter gamma from (1) or (2)
%  - _:     The third argument is a vector y_T or a function handle y_d,
%           depending on whether the type is TerminalCost or Tracking,
%           respectively
%
classdef Obj
    properties
        type
        gamma
        y_T % Terminal cost
        y_d % Tracking
    end
    
    methods
        function obj = Obj(type, gamma, varargin)
            obj.type = type;
            obj.gamma = gamma;
            
            assert(numel(varargin) == 1)
            
            switch type
                case ObjType.TerminalCost, obj.y_T = varargin{1};
                case ObjType.Tracking, obj.y_d = varargin{1};
            end
        end
    end
end

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
            switch type
                case ObjType.TerminalCost, obj.y_T = varargin{1};
                case ObjType.Tracking, obj.y_d = varargin{1};
            end
        end
    end
end

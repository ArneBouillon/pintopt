classdef Obj
    properties
        type
        gamma
        y_d
        y_track
    end
    
    methods
        function obj = Obj(type, gamma, varargin)
            obj.type = type;
            obj.gamma = gamma;
            switch type
                case ObjType.Tracking, obj.y_track = varargin{1};
                case ObjType.TerminalCost, obj.y_d = varargin{1};
            end
        end
    end
end

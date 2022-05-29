classdef Prec
    properties
        type
        alpha
        test
    end

    methods
        function prec = Prec(type, props, obj)
            prec.type = type;
            
            assert(isempty(setdiff(fieldnames(props), {'alpha' 'test'})))
            if isfield(props, 'alpha'), prec.alpha = props.alpha; else, prec.alpha = 1; end
            if isfield(props, 'test'), prec.test = props.test; else, prec.test = false; end

            switch type
                case PrecType.Square
                    assert(abs(abs(prec.alpha)-1) < 1e-8)
                case PrecType.Triangular
                    assert(obj.type == ObjType.TerminalCost)
            end
        end
    end
end

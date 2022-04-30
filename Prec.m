classdef Prec
    properties
        type
        alpha
        diag
        test
    end

    methods
        function prec = Prec(type, props, obj)
            prec.type = type;

            switch type
                case PrecType.Tracking
                    assert(obj.type == ObjType.Tracking)
                    assert(isempty(setdiff(fieldnames(props), {'alpha' 'diag' 'test'})))
                    if isfield(props, 'alpha'), prec.alpha = props.alpha; else, prec.alpha = 1; end
                    if isfield(props, 'test'), prec.test = props.test; else, prec.test = false; end
                    if isfield(props, 'diag'), prec.diag = props.diag; else, prec.diag = false; end
                case PrecType.TerminalCost1
                    assert(obj.type == ObjType.TerminalCost)
                    assert(isempty(setdiff(fieldnames(props), {'alpha' 'test'})))
                    if isfield(props, 'alpha'), prec.alpha = props.alpha; else, prec.alpha = 1; end
                    if isfield(props, 'test'), prec.test = props.test; else, prec.test = false; end
                    prec.diag = false;
                case PrecType.TerminalCost2
                    assert(obj.type == ObjType.TerminalCost)
                    assert(isempty(setdiff(fieldnames(props), {'alpha' 'test'})))
                    if isfield(props, 'alpha'), prec.alpha = props.alpha; else, prec.alpha = 1; end
                    if isfield(props, 'test'), prec.test = props.test; else, prec.test = false; end
                    prec.diag = false;
                case PrecType.TerminalCostMod
                    assert(obj.type == ObjType.TerminalCostMod)
                    assert(isempty(setdiff(fieldnames(props), {'alpha' 'diag' 'test'})))
                    if isfield(props, 'alpha'), prec.alpha = props.alpha; else, prec.alpha = 1; end
                    if isfield(props, 'test'), prec.test = props.test; else, prec.test = false; end
                    if isfield(props, 'diag'), prec.diag = props.diag; else, prec.diag = false; end
            end
        end
    end
end

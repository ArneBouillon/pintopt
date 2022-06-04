%% Preconditioner
% This class contains information about ParaDiag or ParaOpt
% preconditioners. It has the following properties.
%  - type:  The type of preconditioner; see PrecType.m
%  - props: A structure containing properties of the preconditioner.
%           Currently allowed are
%            - alpha: The alpha value of the block-alpha-circulant
%                     preconditioner
%            - test:  When set to true, the algorithm used is obligated to
%                     implement the precontioner as a (sparse) matrix,
%                     instead of a procedure that implicitly solves a
%                     system. This can be used in debugging.
%  - obj:  The Obj instance used for the problem -- used to perform checks
%
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

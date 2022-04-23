classdef Prec
    properties
        alpha
        diag
        test
    end

    methods
        function prec = Prec(alpha, diag, test)
            if ~exist('test', 'var') || isempty(test), test = false; end

            prec.alpha = alpha;
            prec.diag = diag;
            prec.test = test;
        end
    end
end

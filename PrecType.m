%% Preconditioner types
% An enumeration of preconditioner types, containing
%  - Square:     Representing the preconditioner as a 2x2 block matrix, all
%                four blocks may be non-zero and inverting the preconditioner
%                couples forward and backward equations. Alpha must have a
%                magnitude of 1.
%  - Triangular: In the same representation as above, the bottom-left block
%                is emptied to form a preconditioner, such that the
%                backward part can be solved first.
%
%  For more info on the preconditioner types, consult Chapters 4 and 5 of
%  the thesis.
%
classdef PrecType
    enumeration
        Square, Triangular
    end
end

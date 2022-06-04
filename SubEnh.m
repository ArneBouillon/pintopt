%% Subspace enhancement
% Enumeration which contains the various types of subspace enhancement that
% can be applied to ParaOpt (see paraopt.m). The options are
%  - None:        No subspace enhancement at all
%  - Generic:     Subspace enhancement which considers the forward and backward
%                 ParaOpt propagators seperately (NOTE: ensure the fine
%                 propagators are affine in their state and adjoint arguments
%                 before using this type)
%  - Specialized: Subspace enhancement that takes the typical relationship
%                 between forward and backward propagators into account.
%                 See §6.6.2 of the thesis for more info. The NOTE from
%                 generic subspace enhancement also applies here.
classdef SubEnh
    properties
        any
    end
    
   methods
      function subenh = SubEnh(any)
         subenh.any = any;
      end
   end
    
    enumeration
        None(false), Generic(true), Specialized(true)
    end
end

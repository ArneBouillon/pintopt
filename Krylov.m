classdef Krylov
    properties
        any
    end
    
   methods
      function krylov = Krylov(any)
         krylov.any = any;
      end
   end
    
    enumeration
        None(false), Generic(true), Specialized(true)
    end
end

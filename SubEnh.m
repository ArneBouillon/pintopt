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

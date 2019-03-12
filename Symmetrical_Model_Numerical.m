
classdef Symmetrical_Model_Numerical(xu, xa, x0, zu, za, z0, zq, mu, ma, m0, mq)
   properties
      A = [[xu, xa, x0, 0],[zu, za, z0, zq], [0, 0, 0, V/c], [mu, ma, m0, mq]]
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value] * n;
      end
   end
end


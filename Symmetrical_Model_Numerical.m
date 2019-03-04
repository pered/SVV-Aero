
classdef Symmetrical_Model_Numerical
   properties
      A = [[xu, xalpha, xtheta, 0],[zu, zalpha, ztheta, zq], [0, 0, 0, V/c], [mu, malpha, mtheta, mq]]
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
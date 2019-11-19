function tests = ValveTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)  
cd Tests
end

function TestLinear(testCase)
Valve_Test_Linear;
end

function TestNonlinear(testCase)
Valve_Test_Nonlinear;
end

function TestDerivativeSA(testCase)
Valve_Test_Derivative_dSA;
end

function TestDerivativeJ(testCase)
Valve_Test_Derivative_dJ;
end




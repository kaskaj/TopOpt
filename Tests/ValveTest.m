function tests = ValveTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    close all;
    cd ..
end

function TestLinear(testCase)
    Valve_Test_Linear;
end

function TestNonlinear(testCase)
    Valve_Test_Nonlinear;
end

function TestInitial(testCase)
    Valve_Test_Initial;
end

function TestDerivativeSA(testCase)
    Valve_Test_Derivative_dJ_nonlinear;
end

function TestDerivativeJ(testCase)
    Valve_Test_Derivative_dJ;
end




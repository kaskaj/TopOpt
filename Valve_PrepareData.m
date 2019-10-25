levels_max = 3;

% Set parameters
Valve_SetParameters;

% Create topology
Valve_MeshGen();

% Create matrices
for levels = 1:levels_max
    fprintf('Creating matrix level %d.\n', levels);
    if levels <= 3
        Valve_CreateMatrices(levels,1);
    else
        Valve_CreateMatrices(levels,1);
    end
end


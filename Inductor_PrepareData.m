levels_max = 5;

% Set parameters
Inductor_SetParameters;

% Create topology
Inductor_MeshGen();

% Create matrices
for levels = 1:levels_max
    fprintf('Creating matrix level %d.\n', levels);
    if levels <= 3
        Inductor_CreateMatrices(levels,1);
    else
        Inductor_CreateMatrices(levels,1);
    end
end


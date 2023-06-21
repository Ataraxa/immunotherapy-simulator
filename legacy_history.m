% Initial conditions for the model
function s = legacy_history(t)
    % Need to to ask Dr. Tanaka what are baseline biological levels
    s = [
        0.025; % Amount of IFNg
        12; % Amount of CD8+
        4.1; % Amount of PD-1
        7.12; % Volume of injected living tumour
        0; % Initial volume of dead tumour
    ];
end
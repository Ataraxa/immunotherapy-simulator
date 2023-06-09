% Initial conditions for the model
function s = history(t)
    s = [
        0; % Amount of IFNg
        0; % Amount of CD8+
        0; % Amount of PD-1
        7.12; % Volume of injected living tumour
        0; % Initial volume of dead tumour
    ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for comlpex refractive indices look in: https://refractiveindex.info/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function properties = get_material_properties(material_type)

switch material_type
    case 'gold'
        alphaT = 14.2e-6; % linear thermal expansion coefficient of sample (m/(m*degC))
        C = 130; % specific heat capacity of sample (J/(kg*degC))
        rho = 19290; % mass density of sample (kg/m^3)
        kappa = 315; % thermal conductivity of sample (W/(m*degC))
        E = 74e9; % Young's modulus (Pa)
        n_520nm = 0.57669+2.1781i; % complex refractive index at 520 nm
        
    case 'silicon'
        alphaT = 4e-6; % linear thermal expansion coefficient of sample (m/(m*degC))
        C = 710; % specific heat capacity of sample (J/(kg*degC))
        rho = 2330; % mass density of sample (kg/m^3)
        kappa = 154; % thermal conductivity of sample (W/(m*degC))
        E = 110e9; % Young's modulus (Pa)
        n_520nm = 4.1870+0.039531i; % complex refractive index at 520 nm
        
    case 'organic dye'
        alphaT = 70e-6; % linear thermal expansion coefficient of sample (m/(m*degC))
        C = 1670; % specific heat capacity of sample (J/(kg*degC))
        rho = 1030; % mass density of sample (kg/m^3)
        kappa = 0.03; % thermal conductivity of sample (W/(m*degC))
        E = 3e9; % Young's modulus (Pa)
        n_520nm = 3.5+2i; % complex refractive index at 520 nm
        
    case 'glass'
        alphaT = 6e-6; % linear thermal expansion coefficient of sample (m/(m*degC))
        C = 840; % specific heat capacity of sample (J/(kg*degC))
        rho = 2500; % mass density of sample (kg/m^3)
        kappa = 1.05; % thermal conductivity of sample (W/(m*degC))
        E = 70e9; % Young's modulus (Pa)
        n_520nm = 1.5202+8.4423e-9i;  % complex refractive index at 520 nm
        
    otherwise
        alphaT = 6e-6; % linear thermal expansion coefficient of sample (m/(m*degC))
        C = 840; % specific heat capacity of sample (J/(kg*degC))
        rho = 2500; % mass density of sample (kg/m^3)
        kappa = 1.05; % thermal conductivity of sample (W/(m*degC))
        E = 70e9; % Young's modulus (Pa)
        n_520nm = 1.5202+8.4423e-9i;  % complex refractive index at 520 nm
        
end
properties = [alphaT; C; rho; kappa; E; n_520nm];

end


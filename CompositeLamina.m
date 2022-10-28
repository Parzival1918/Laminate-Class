classdef CompositeLamina < handle
    %COMPOSITELAMINA Summary: Use this class to create a lamina of aligned
    %fibre reinforcement in a matrix.
    %   Detailed explanation goes here
    
    properties (SetAccess = private)
        S (6,6) {mustBeFloat} = zeros(6)
        C (6,6) {mustBeFloat} = zeros(6)

        laminaE1 (1,1) double {mustBeFloat} = 0
        laminaE2 (1,1) double {mustBeFloat} = 0
        laminav12 (1,1) double {mustBeFloat} = 0
        laminav21 (1,1) double {mustBeFloat} = 0
        laminav23 (1,1) double {mustBeFloat} = 0
        laminaG12 (1,1) double {mustBeFloat} = 0
        laminaG13 (1,1) double {mustBeFloat} = 0
        laminaK (1,1) double {mustBeFloat} = 0
        laminaIR121 (1,1) double {mustBeFloat} = 0
        laminaIR122 (1,1) double {mustBeFloat} = 0
    end

    properties
        compositeComponents CompositeComponents = CompositeComponents('E-glass fibres');
        fibre_fraction (1,1) double {mustBeFloat} = 0.5
        Halpin_Tsai_EqnParameter (1,1) double {mustBeFloat} = 1

        lamina_type char {mustBeMember(lamina_type,{'isotropic',...
                    'orthotropic','transversely isotropic'})} = 'transversely isotropic'
        loading_type char {mustBeMember(loading_type,{'plane stress',...
                    })} = 'plane stress'
    end
    
    methods
        function obj = CompositeLamina(CC)
            %COMPOSITELAMINA Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                CC CompositeComponents = CompositeComponents('E-glass fibres');
            end

            obj.compositeComponents = CC;
        end

        function laminaE1 = get.laminaE1(obj)
            Ef = obj.compositeComponents.reinforcementE1;
            Em = obj.compositeComponents.matrixE1;
            ff = obj.fibre_fraction;

            laminaE1 = Ef*ff + (1-ff)*Em;
        end

        function laminaE2 = get.laminaE2(obj)
            Ef = obj.compositeComponents.reinforcementE1;
            Em = obj.compositeComponents.matrixE1;
            e = obj.Halpin_Tsai_EqnParameter;
            nu = (Ef-Em)/(Ef+e*Em);
            ff = obj.fibre_fraction;

            laminaE2 = Em*((1+e*nu*ff)/(1-nnu*ff));
        end

        function laminaG12 = get.laminaG12(obj)
            Gf = obj.compositeComponents.reinforcementG12;
            Gm = obj.compositeComponents.matrixG12;
            e = obj.Halpin_Tsai_EqnParameter;
            nu = (Gf-Gm)/(Gf+e*Gm);
            ff = obj.fibre_fraction;

            laminaG12 = Gm*((1+e*nu*ff)/(1-nnu*ff));
        end

        function laminav12 = get.laminav12(obj)
            vf = obj.compositeComponents.reinforcementv12;
            vm = obj.compositeComponents.matrixv12;
            ff = obj.fibre_fraction;

            laminav12 = vf*ff + (1-ff)*vm;
        end

        function laminav21 = get.laminav21(obj)
            laminav21 = obj.laminaE2*(obj.laminav12/obj.laminaE1);
        end

        function laminaK = get.laminaK(obj)
            Kf = obj.compositeComponents.reinforcementK;
            Km = obj.compositeComponents.matrixK;
            ff = obj.fibre_fraction;

            laminaK = 1/(ff/Kf + (1-ff)/Km);
        end

        function laminav23 = get.laminav23(obj)
            v21 = obj.laminav21;
            E2 = obj.laminaE2;
            K = obj.laminaK;

            laminav23 = 1 - v21 - E2/(3*K);
        end

        function S = get.S(obj)
            switch obj.lamina_type
                case 'isotropic'
                    E1 = obj.laminaE1;
                    v12 = obj.laminav12;
                    S11 = 1/E1;
                    S12 = -v12/E1;
                    S11S12 = 2*(S11-S12);
                    
                    S = [S11 S12 S12 0 0 0; ...
                         S12 S11 S12 0 0 0; ...
                         S12 S12 S11 0 0 0; ...
                         0 0 0 S11S12 0 0; ...
                         0 0 0 0 S11S12 0; ...
                         0 0 0 0 0 S11S12];
                case 'orthotropic'
                    S = zeros(6);
                case 'transversely isotropic'
                    E1 = obj.laminaE1;
                    E2 = obj.laminaE2;
                    v12 = obj.laminav12;
                    G12 = obj.laminaG12;
                    S11 = 1/E1;
                    S22 = 1/E2;
                    S12 = -v12/E1;
                    S66 = 1/G12;

                    S = [S11 S12 0 0 0 0; ...
                         S12 S22 0 0 0 0; ...
                         0 0 0 0 0 0; ...
                         0 0 0 0 0 0; ...
                         0 0 0 0 0 0; ...
                         0 0 0 0 0 S66];
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

function val = U(angle)
    val = [cosd(angle)^2 sind(angle)^2 cosd(angle)*sind(angle); ...
           sind(angle)^2 cosd(angle)^2 -cosd(angle)*sind(angle); ...
           -2*cosd(angle)*sind(angle) 2*cosd(angle)*sind(angle) cosd(angle)^2 - sind(angle)^2];
end
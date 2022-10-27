classdef Lamina < Composite
    properties (SetAccess = private)
        stress (6,1) {mustBeFloat} = zeros([6 1])
        strain (6,1) {mustBeFloat} = zeros([6 1])
        C (6,6) {mustBeFloat} = zeros(6)
        S (6,6) {mustBeFloat} = zeros(6)
        material char {mustBeMember(material,{'isotropic', 'orthotropic,transversely isotropic'})} ...
            = 'orthotropic,transversely isotropic'
        angle_Z (1,1) {mustBeFloat} = 0
    end
    
    methods

        function obj = Lamina(Efibres, Ematrix, vFibre, vMatrix)
            if nargin == 1
                obj.FibreE = Efibres;
            elseif nargin == 2
                obj.FibreE = Efibres;
                obj.MatrixE = Ematrix;
            elseif nargin == 3
                obj.FibreE = Efibres;
                obj.MatrixE = Ematrix;
                obj.Fibrev = vFibre;
            elseif nargin == 4
                obj.FibreE = Efibres;
                obj.MatrixE = Ematrix;
                obj.Fibrev = vFibre;
                obj.Matrixv = vMatrix;
                obj.update_Composite;
            end
        end

        function obj = new_stress_material(obj,stress)
            arguments
                obj {nonZeroS(obj)}
                stress {mustBeOfMaterialType(obj,stress)}
            end

            obj.stress = stress;
            obj.update_strain(obj.S,obj.stress);
        end

        function obj = new_stress_angledZ(obj,stress,angleZ)
            arguments
                obj {nonZeroS(obj)}
                stress {mustBeOfMaterialType(obj,stress)}
                angleZ {mustBeFloat}
            end

            if size(stress,1) == 6
                stress = to_tensor(stress);
            elseif ne(size(stress,1),3) && ne(size(stress,2),3)
                error("Input stress tensor of wrong dimensions.")
            end

            equals = obj.S == zeros(6);
        
            if equals == numel(zeros(6))
                error("Set a valid value for compliance matrix.")
            end

            %rotation_matrix = [cosd(angleZ) -sind(angleZ) 0; ...
            %                   sind(angleZ) cosd(angleZ) 0; ...
            %                   0 0 1];
            cossin = [cosd(angleZ)^2 sind(angleZ)^2 cosd(angleZ)*sind(angleZ); ...
                          sind(angleZ)^2 cosd(angleZ)^2 -cosd(angleZ)*sind(angleZ); ...
                          -2*cosd(angleZ)*sind(angleZ) 2*cosd(angleZ)*sind(angleZ) cosd(angleZ)^2 - sind(angleZ)^2];
                        
            rotated_stress = cossin*stress;
            obj.stress = to_compact(rotated_stress);

            obj.update_strain(obj.S,obj.stress);
            obj.angle_Z = angleZ;
        end

        function obj = update_Stiffness_Compliance(obj)
            if obj.FibreE == 0 || obj.FibreG == 0 || obj.Fibrev == 0 || ...
                    obj.MatrixE == 0 || obj.MatrixG == 0 || obj.Matrixv == 0
                error("Composite properties of fibre and matrix are not all updated.")
            end
            
            if obj.material == "isotropic"
                obj.new_S(obj.compositeE1,obj.compositev12);
            elseif obj.material == "orthotropic,transversely isotropic"
                obj.new_S(obj.compositeE1,obj.compositev12,obj.compositeE2,obj.compositeG12);
            end
        end

        function obj = set_material(obj,prop)
            obj.material = prop;

            obj.S = zeros(6);
            obj.C = zeros(6);
        end

        function YMs = composite_YMs(obj)
            if obj.material == "orthotropic,transversely isotropic"
                E1 = 1/obj.S(1,1);
                E2 = 1/obj.S(2,2);
                YMs = [E1 E2];
            else
                E1 = 1/obj.S(1,1);
                YMs = E1;
            end
        end

        function YMs = angle_YMs(obj)
            if obj.material == "orthotropic,transversely isotropic"
                angleZ = obj.angle_Z;
                cossin = [cosd(angleZ)^2 sind(angleZ)^2 cosd(angleZ)*sind(angleZ); ...
                          sind(angleZ)^2 cosd(angleZ)^2 -cosd(angleZ)*sind(angleZ); ...
                          -2*cosd(angleZ)*sind(angleZ) 2*cosd(angleZ)*sind(angleZ) cosd(angleZ)^2 - sind(angleZ)^2];
                        
                S_ = [obj.S(1,1), obj.S(1,2), 0; ...
                     obj.S(2,1), obj.S(2,2), 0; ...
                     0, 0, obj.S(6,6)];
                S_ = cossin*S_*transpose(cossin);
                E1 = 1/S_(1,1);
                E2 = 1/S_(2,2);
                YMs = [E1 E2];
            else
                E1 = 1/obj.S(1,1);
                YMs = E1;
            end
        end

        function G = composite_G(obj)
            if obj.material == "orthotropic,transversely isotropic"
                G12 = 1/obj.S(6,6);
                G = G12;
            else
                G = NaN;
            end
        end

        function G = angle_G(obj)
            if obj.material == "orthotropic,transversely isotropic"
                angleZ = obj.angle_Z;
                cossin = [cosd(angleZ)^2 sind(angleZ)^2 cosd(angleZ)*sind(angleZ); ...
                          sind(angleZ)^2 cosd(angleZ)^2 -cosd(angleZ)*sind(angleZ); ...
                          -2*cosd(angleZ)*sind(angleZ) 2*cosd(angleZ)*sind(angleZ) cosd(angleZ)^2 - sind(angleZ)^2];
                        
                S_ = [obj.S(1,1), obj.S(1,2), 0; ...
                     obj.S(2,1), obj.S(2,2), 0; ...
                     0, 0, obj.S(6,6)];
                S_ = cossin*S_*transpose(cossin);
                G12 = 1/S_(3,3);
                G = G12;
            else
                G = NaN;
            end
        end

        function Vs = composite_v(obj)
            if obj.material == "orthotropic,transversely isotropic"
                E1 = 1/obj.S(1,1);
                E2 = 1/obj.S(2,2);
                v12 = -E1*obj.S(1,2);
                v21 = -E2*obj.S(1,2);
                Vs = [v12 v21];
            else
                E1 = 1/obj.S(1,1);
                v12 = -E1*obj.S(1,2);
                Vs = v12;
            end
        end

        function Vs = angle_v(obj)
            if obj.material == "orthotropic,transversely isotropic"
                angleZ = obj.angle_Z;
                cossin = [cosd(angleZ)^2 sind(angleZ)^2 cosd(angleZ)*sind(angleZ); ...
                          sind(angleZ)^2 cosd(angleZ)^2 -cosd(angleZ)*sind(angleZ); ...
                          -2*cosd(angleZ)*sind(angleZ) 2*cosd(angleZ)*sind(angleZ) cosd(angleZ)^2 - sind(angleZ)^2];
                S_ = [obj.S(1,1), obj.S(1,2), 0; ...
                     obj.S(2,1), obj.S(2,2), 0; ...
                     0, 0, obj.S(6,6)];
                S_ = cossin*S_*transpose(cossin);

                E1 = 1/S_(1,1);
                E2 = 1/S_(2,2);
                v12 = -E1*S_(1,2);
                v21 = -E2*S_(1,2);
                Vs = [v12 v21];
            else
                E1 = 1/obj.S(1,1);
                v12 = -E1*obj.S(1,2);
                Vs = v12;
            end
        end

        function interaction_ratios = composite_interaction_ratios(obj)
            if obj.material == "orthotropic,transversely isotropic"
                angleZ = obj.angle_Z;
                E1 = 1/obj.S(1,1);
                E2 = 1/obj.S(2,2);
                v12 = -E1*obj.S(1,2);
                G12 = 1/obj.S(6,6);

                S16 = (sind(2*angleZ)*(cosd(angleZ)^2))/E1 - (sind(2*angleZ)*(sind(angleZ)^2))/E2 + ...
                      (sind(4*angleZ)*v12)/(2*E1) - sind(4*angleZ)/(4*G12);

                S26 = (sind(2*angleZ)*(sind(angleZ)^2))/E1 - (sind(2*angleZ)*(cosd(angleZ)^2))/E2 - ...
                      (sind(4*angleZ)*v12)/(2*E1) + sind(4*angleZ)/(4*G12);

                ratio121 = S16*E1;

                ratio122 = S26*E2;

                interaction_ratios = [ratio121 ratio122];
            else
                interaction_ratios = NaN;
            end
        end

        function obj = new_angleZ(obj, angle)
            obj.angle_Z = angle;
        end

        function plot_YMs_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end

            YMsAngle = zeros([length(angles), 2]);
            pos = 1;
            for angle = angles
                obj.new_angleZ(angle);
                YMsAngle(pos,:) = obj.angle_YMs();
                pos = pos + 1;
            end
            
            figure;
            plot(angles,YMsAngle);
            xlabel("Angle (º)");
            ylabel("Major and Minor Modulus (" + obj.stiffness_units + ")");
            legend("E1","E2");
        end

        function plot_G_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end

            GAngle = zeros([length(angles), 1]);
            pos = 1;
            for angle = angles
                obj.new_angleZ(angle);
                GAngle(pos,:) = obj.angle_G();
                pos = pos + 1;
            end
            
            figure;
            plot(angles,GAngle);
            xlabel("Angle (º)");
            ylabel("G12 (" + obj.stiffness_units + ")");
        end

        function plot_v_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end

            vAngle = zeros([length(angles), 2]);
            pos = 1;
            for angle = angles
                obj.new_angleZ(angle);
                vAngle(pos,:) = obj.angle_v();
                pos = pos + 1;
            end
            
            figure;
            plot(angles,vAngle);
            xlabel("Angle (º)");
            ylabel("Major and Minor Poisson's Ratio");
            legend("v12","v21");
        end

        function plot_interaction_ratios_angle(obj,angles)
            arguments
                obj
                angles (1,:) {mustBeFloat} = linspace(0,90,30)
            end

            IRAngle = zeros([length(angles), 2]);
            pos = 1;
            for angle = angles
                obj.new_angleZ(angle);
                IRAngle(pos,:) = obj.composite_interaction_ratios();
                pos = pos + 1;
            end
            
            figure;
            plot(angles,IRAngle);
            xlabel("Angle (º)");
            ylabel("Interaction Ratios");
            legend("n121","n122");
        end
    end
    methods (Access = private)
        function obj = update_strain(obj,S,stress)
            obj.strain = S*(stress*obj.ratio);
        end

        function obj = new_S(obj,varargin)
            if obj.material == "isotropic" && nargin == 3
                S11 = 1/varargin{1};
                S12 = -varargin{2}/varargin{1};
                obj.S = [S11, S12, S12, 0, 0, 0; ...
                         S12, S11, S12, 0, 0, 0; ...
                         S12, S12, S11, 0, 0, 0; ...
                         0, 0, 0, 2*(S11-S12), 0, 0; ...
                         0, 0, 0, 0, 2*(S11-S12), 0; ...
                         0, 0, 0, 0, 0, 2*(S11-S12)];
                obj.C = inv(obj.S);
            elseif obj.material == "orthotropic,transversely isotropic" && nargin == 5
                S11 = 1/varargin{1};
                S12 = -varargin{2}/varargin{1};
                S22 = 1/varargin{3};
                S66 = 1/varargin{4};
                obj.S = [S11, S12, 0, 0, 0, 0; ...
                         S12, S22, 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, S66];
                tempS = [S11, S12, 0; ...
                         S12, S22, 0; ...
                         0, 0, S66];
                tempC = inv(tempS);
                obj.C = [tempC(1,1), tempC(1,2), 0, 0, 0, 0; ...
                         tempC(2,1), tempC(2,2), 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, 0; ...
                         0, 0, 0, 0, 0, tempC(3,3)]; %%ERROR
            else
                error("invalid number of inputs for Laminate material type.")
            end
        end
    end
end

function mustBeOfMaterialType(obj, stress)
    if size(stress,1) == 3 && size(stress,2) == 3
        stress = to_compact(stress);
    end

    material_type = obj.material;
    if material_type == "orthotropic,transversely isotropic" && ...
      (ne(stress(3),0) || ne(stress(4),0) || ne(stress(5),0))
        error("Invalid applied stress for material properties assigned. positions 3, 4, 5 must be 0.")
    end
end

function nonZeroS(obj)
    equals = obj.S == zeros(6);

    if equals == numel(zeros(6))
        error("Set a valid value for compliance matrix.")
    end
end

function compact = to_compact(tensor)
    compact = [tensor(1,1); tensor(2,2); tensor(3,3); tensor(2,3); ...
               tensor(1,3); tensor(1,2)];
end

function tensor = to_tensor(compact)
    tensor = [compact(1) compact(6) compact(5); ...
              compact(6) compact(2) compact(4); ...
              compact(5) compact(4) compact(3)];
end
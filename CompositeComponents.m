classdef CompositeComponents < handle
    %COMPOSITECOMPONENTS Summary: Use this class to add material properties
    %of the fibre and matrix materials. Then these properties can be used
    %to create a lamina.
    %   Detailed explanation goes here

    properties 
        reinforcementE1 (1,1) double {mustBeFloat} = 0 %Reinforcement young's 
                                                % modulus
        matrixE1 (1,1) double {mustBeFloat} = 0
        reinforcementv12 (1,1) double {mustBeFloat} = 0
        matrixv12 (1,1) double {mustBeFloat} = 0
        reinforcementG12 (1,1) double {mustBeFloat} = 0
        matrixG12 (1,1) double {mustBeFloat} = 0
        reinforcementK (1,1) double {mustBeFloat} = 0 %Reinforcement bulk modulus
        matrixK (1,1) double {mustBeFloat} = 0
        
        reinforcement_density (1,1) double {mustBeFloat} = 0
        matrix_denisty (1,1) double {mustBeFloat} = 0

        E_units char {mustBeMember(E_units,{'GPa','MPa','Pa'})} = 'GPa'
        v_units char {mustBeMember(v_units,{'unitless'})} = 'unitless'
        G_units char {mustBeMember(G_units,{'GPa','MPa','Pa'})} = 'GPa'
        K_units char {mustBeMember(K_units,{'GPa','MPa','Pa'})} = 'GPa'
        density_units char {mustBeMember(density_units,{'g/cm3'})} = 'g/cm3'

        reinforcement_type char {mustBeMember(reinforcement_type, ...
                           {'fibre'})} = 'fibre'
        matrix_type char {mustBeMember(matrix_type, ...
                    {'thermoset','thermoplastic'})} = 'thermoset'
    end
    
    methods
        function obj = CompositeComponents(varargin)
            %COMPOSITECOMPONENTS Construct an instance of this class
            %   Detailed explanation goes here
            narginchk(0,8);
            switch nargin
                case 1
                    if isa(varargin{1},'double')
                        obj.reinforcementE1 = varargin{1};
                    elseif ischar(class(varargin{1}))
                        obj.assign_reinforcement_material(varargin{1})
                    end
                case 2
                    obj.reinforcementE1 = varargin{1};
                    obj.matrixE1 = varargin{2};
                case 3
                    obj.reinforcementE1 = varargin{1};
                    obj.matrixE1 = varargin{2};
                    obj.reinforcementv12 = varargin{3};
                case 4
                    obj.reinforcementE1 = varargin{1};
                    obj.matrixE1 = varargin{2};
                    obj.reinforcementv12 = varargin{3};
                    obj.matrixv12 = varargin{4};
                case 5
                    obj.reinforcementE1 = varargin{1};
                    obj.matrixE1 = varargin{2};
                    obj.reinforcementv12 = varargin{3};
                    obj.matrixv12 = varargin{4};
                    obj.reinforcementG12 = varargin{5};
                case 6
                    obj.reinforcementE1 = varargin{1};
                    obj.matrixE1 = varargin{2};
                    obj.reinforcementv12 = varargin{3};
                    obj.matrixv12 = varargin{4};
                    obj.reinforcementG12 = varargin{5};
                    obj.matrixG12 = varargin{6};
                case 7
                    obj.reinforcementE1 = varargin{1};
                    obj.matrixE1 = varargin{2};
                    obj.reinforcementv12 = varargin{3};
                    obj.matrixv12 = varargin{4};
                    obj.reinforcementG12 = varargin{5};
                    obj.matrixG12 = varargin{6};
                    obj.reinforcementK = varargin{7};
                case 8
                    obj.reinforcementE1 = varargin{1};
                    obj.matrixE1 = varargin{2};
                    obj.reinforcementv12 = varargin{3};
                    obj.matrixv12 = varargin{4};
                    obj.reinforcementG12 = varargin{5};
                    obj.matrixG12 = varargin{6};
                    obj.reinforcementK = varargin{7};
                    obj.matrixK = varargin{8};
            end
        end
        
        function obj = calc_property(obj, fibre_or_matrix, property)
            %METHOD1 Summary: Update one of the material properties using
            %values of other ones
            %   Detailed explanation goes here
            arguments
                obj CompositeComponents
                fibre_or_matrix {mustBeMember(fibre_or_matrix, ...
                                {'fibre','matrix'})} = 'fibre'
                property {mustBeMember(property,{'E','v','G','K'})} = 'G'
            end
            
            if strcmp(fibre_or_matrix,'fibre')
                if strcmp(property,'E')
                    if ne(obj.reinforcementG12,0) && ...
                       ne(obj.reinforcementv12,0)
                        obj.reinforcementE1 = 2*obj.reinforcementG12* ...
                                              (1 + obj.reinforcementv12);
                    elseif ne(obj.reinforcementK,0) && ...
                           ne(obj.reinforcementv12,0)
                        obj.reinforcementE1 = 3*obj.reinforcementK* ...
                                              (1 - 2*obj.reinforcementv12);
                    end
                elseif strcmp(property,'v')
                    if ne(obj.reinforcementG12,0) && ...
                       ne(obj.reinforcementE1,0)
                        obj.reinforcementv12 = obj.reinforcementE1/ ...
                                               (obj.reinforcementG12*2) -1;
                    elseif ne(obj.reinforcementK,0) && ...
                           ne(obj.reinforcementE1,0)
                        obj.reinforcementv12 = 0.5 - obj.reinforcementE1/ ...
                                               (6*obj.reinforcementK);
                    end
                elseif strcmp(property,'G')
                    if ne(obj.reinforcementv12,0) && ...
                       ne(obj.reinforcementE1,0)
                        obj.reinforcementG12 = obj.reinforcementE1/ ...
                                               (2*(1 + obj.reinforcementv12));
                    elseif ne(obj.reinforcementv12,0) && ...
                           ne(obj.reinforcementK,0)
                        obj.reinforcementG12 = (3*obj.reinforcementK(1-...
                                               2*obj.reinforcementv12))/ ...
                                               (2*(1 + obj.reinforcementv12));
                    end
                else
                    if ne(obj.reinforcementv12,0) && ...
                       ne(obj.reinforcementE1,0)
                        obj.reinforcementK = obj.reinforcementE1/ ...
                                             (3*(1-2*obj.reinforcementv12));
                    elseif ne(obj.reinforcementv12,0) && ...
                           ne(obj.reinforcementG12,0)
                        obj.reinforcementK = (obj.reinforcementG12* ...
                                             (2*(1 + obj.reinforcementv12)))/ ...
                                             (3*(1-2*obj.reinforcementv12));
                    end
                end
            else
                if strcmp(property,'E')
                    if ne(obj.matrixG12,0) && ...
                       ne(obj.matrixv12,0)
                        obj.matrixE1 = 2*obj.matrixG12* ...
                                              (1 + obj.matrixv12);
                    elseif ne(obj.matrixK,0) && ...
                           ne(obj.matrixv12,0)
                        obj.matrixE1 = 3*obj.matrixK* ...
                                              (1 - 2*obj.matrixv12);
                    end
                elseif strcmp(property,'v')
                    if ne(obj.matrixG12,0) && ...
                       ne(obj.matrixE1,0)
                        obj.matrixv12 = obj.matrixE1/ ...
                                               (obj.matrixG12*2) -1;
                    elseif ne(obj.matrixK,0) && ...
                           ne(obj.matrixE1,0)
                        obj.matrixv12 = 0.5 - obj.matrixE1/ ...
                                               (6*obj.matrixK);
                    end
                elseif strcmp(property,'G')
                    if ne(obj.matrixv12,0) && ...
                       ne(obj.matrixE1,0)
                        obj.matrixG12 = obj.matrixE1/ ...
                                               (2*(1 + obj.matrixv12));
                    elseif ne(obj.matrixv12,0) && ...
                           ne(obj.matrixK,0)
                        obj.matrixG12 = (3*obj.matrixK(1-...
                                               2*obj.matrixv12))/ ...
                                               (2*(1 + obj.matrixv12));
                    end
                else
                    if ne(obj.matrixv12,0) && ...
                       ne(obj.matrixE1,0)
                        obj.matrixK = obj.matrixE1/ ...
                                             (3*(1-2*obj.matrixv12));
                    elseif ne(obj.matrixv12,0) && ...
                           ne(obj.matrixG12,0)
                        obj.matrixK = (obj.matrixG12* ...
                                             (2*(1 + obj.matrixv12)))/ ...
                                             (3*(1-2*obj.matrixv12));
                    end
                end
            end
        end

        function obj = assign_reinforcement_material(obj, name)
            %METHOD1 Summary: Assign material properties of the
            %reinforcement from existing materials in the databse
            %   Detailed explanation goes here
            arguments
                obj CompositeComponents
                name {mustBeMember(name,{'carbon fibres','A-glass fibres',...
                     'C-glass fibres','E-glass fibres','S2-glass fibres'})}
            end
            
            switch name
                case 'A-glass fibres'
                    %https://textilelearner.net/glass-fiber-types-properties/
                    obj.reinforcement_type = 'fibre';
                    obj.reinforcementE1 = 72;
                    obj.reinforcement_density = 2.44;
                case 'C-glass fibres'
                    %https://textilelearner.net/glass-fiber-types-properties/
                    obj.reinforcement_type = 'fibre';
                    obj.reinforcementE1 = 69;
                    obj.reinforcement_density = 2.56;
                case 'E-glass fibres'
                    %https://www.researchgate.net/figure/Physical-and-mechanical-properties-of-glass-fiber_tbl2_265346634
                    obj.reinforcement_type = 'fibre';
                    obj.reinforcementE1 = 72.3;
                    obj.reinforcement_density = 2.58;
                    obj.reinforcementv12 = 0.2;
                case 'S2-glass fibres'
                    %https://www.researchgate.net/figure/Physical-and-mechanical-properties-of-glass-fiber_tbl2_265346634
                    obj.reinforcement_type = 'fibre';
                    obj.reinforcementE1 = 86.9;
                    obj.reinforcement_density = 2.46;
                    obj.reinforcementv12 = 0.22;
            end
        end
        
        %{
        function obj = assign_matrix_material(obj, name)
            arguments
                obj CompositeComponents
                name {mustBeMember(name,{})}
            end
            
            switch name
                case 
            end
        end
        %}
    end
end


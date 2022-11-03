classdef CompositeComponents < handle
    %COMPOSITECOMPONENTS Summary: Use this class to add material properties
    %of the fibre and matrix materials. Then these properties can be used
    %to create a lamina.
    %   For now this class is very simple and can only work with materials
    %   that are isotropic. I plan to add support for anisotropic materials
    %   later on.
    
    properties 
        reinforcementE1 (1,1) double {mustBeFloat} = 0 %Reinforcement young's modulus
        matrixE1 (1,1) double {mustBeFloat} = 0 %Matrix young's modulus
        reinforcementv12 (1,1) double {mustBeFloat} = 0 %Reinforcement Poisson's ratio modulus
        matrixv12 (1,1) double {mustBeFloat} = 0 %Matrix Poisson's ratio modulus
        reinforcementG12 (1,1) double {mustBeFloat} = 0 %Reinforcement shear modulus
        matrixG12 (1,1) double {mustBeFloat} = 0 %Matrix shear modulus
        reinforcementK (1,1) double {mustBeFloat} = 0 %Reinforcement bulk modulus
        matrixK (1,1) double {mustBeFloat} = 0 %Matrix bulk modulus
        
        reinforcement_density (1,1) double {mustBeFloat} = 0 %Reinforcement material density
        matrix_denisty (1,1) double {mustBeFloat} = 0 %Matrix material density

        E_units char {mustBeMember(E_units,{'GPa','MPa','Pa'})} = 'GPa' %Units of materials stifness
        v_units char {mustBeMember(v_units,{'unitless'})} = 'unitless' %Units of materials poisson's ratios
        G_units char {mustBeMember(G_units,{'GPa','MPa','Pa'})} = 'GPa' %Units of the shearing moodulus
        K_units char {mustBeMember(K_units,{'GPa','MPa','Pa'})} = 'GPa' %Units of the bulk modulus
        density_units char {mustBeMember(density_units,{'g/cm3'})} = 'g/cm3' %Units of the densities

        reinforcement_type char {mustBeMember(reinforcement_type, ... %Set the type of reinforcement that the composite will ahve
                           {'long fibres','short fibres'})} = 'long fibres' 
        matrix_type char {mustBeMember(matrix_type, ... %Set the type of matrix that the composite will have
                    {'thermoset','thermoplastic'})} = 'thermoset' 
        r (1,1) double {mustBeFloat} = 0.5 %Reinforcement radius
        L (1,1) double {mustBeFloat} = 100 %Reinforcement half-length
    end
    
    methods
        function obj = CompositeComponents(varargin)
            %COMPOSITECOMPONENTS Construct an instance of this class
            %   Can construct the class with some parameters already set,
            %   or call with one char argument out of the accepted ones in 
            %   ASSIGN_REINFORCEMENT_MATERIAL() method. 
            narginchk(0,8);
            switch nargin
                case 1
                    if isa(varargin{1},'double')
                        obj.reinforcementE1 = varargin{1};
                    elseif ischar(class(varargin{1}))
                        obj.assign_reinforcement_material(varargin{1});
                    end
                case 2
                    if isa(varargin{1},'double') && isa(varargin{2},'double')
                        obj.reinforcementE1 = varargin{1};
                        obj.matrixE1 = varargin{2};
                    elseif ischar(class(varargin{1})) && ischar(class(varargin{2}))
                        obj.assign_reinforcement_material(varargin{1},varargin{2});
                    end
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
        
        function obj = calc_property(obj, reinforcement_or_matrix, property)
            %CALC_PROPERTY Summary: Update one of the material properties using
            %values of other ones
            %   Calculate values of materials E, G, v or K using the
            %   formulas (eg E = 2*G*(1 + v)) relating all 4 properties.
            arguments
                obj CompositeComponents
                reinforcement_or_matrix {mustBeMember(reinforcement_or_matrix, ...
                                {'reinforcement','matrix'})} = 'reinforcement'
                property {mustBeMember(property,{'E','v','G','K'})} = 'G'
            end
            
            if strcmp(reinforcement_or_matrix,'reinforcement')
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

        function obj = assign_reinforcement_material(obj, name, length)
            %ASSIGN_REINFORCEMENT_MATERIAL Summary: Assign material properties of the
            %reinforcement from existing materials in the databse
            %   Use the properties of the materials in the database in your
            %   project. I am manually adding them, so maybe an improvement
            %   would be to add some way of scrapping web data or adding
            %   from text files.
            arguments
                obj CompositeComponents
                name {mustBeMember(name,{'carbon fibres','A-glass fibres',...
                     'C-glass fibres','E-glass fibres','S2-glass fibres'})}
                length {mustBeMember(length,{'short fibres','long fibres'})} = 'long fibres'
            end
            
            switch name
                case 'A-glass fibres'
                    %https://textilelearner.net/glass-fiber-types-properties/
                    obj.reinforcement_type = length;
                    obj.reinforcementE1 = 72;
                    obj.reinforcement_density = 2.44;
                case 'C-glass fibres'
                    %https://textilelearner.net/glass-fiber-types-properties/
                    obj.reinforcement_type = length;
                    obj.reinforcementE1 = 69;
                    obj.reinforcement_density = 2.56;
                case 'E-glass fibres'
                    %https://www.researchgate.net/figure/Physical-and-mechanical-properties-of-glass-fiber_tbl2_265346634
                    obj.reinforcement_type = length;
                    obj.reinforcementE1 = 72.3;
                    obj.reinforcement_density = 2.58;
                    obj.reinforcementv12 = 0.2;
                    obj.calc_property('reinforcement','G');
                    obj.calc_property('reinforcement','K');
                case 'S2-glass fibres'
                    %https://www.researchgate.net/figure/Physical-and-mechanical-properties-of-glass-fiber_tbl2_265346634
                    obj.reinforcement_type = length;
                    obj.reinforcementE1 = 86.9;
                    obj.reinforcement_density = 2.46;
                    obj.reinforcementv12 = 0.22;
                    obj.calc_property('reinforcement','G');
                    obj.calc_property('reinforcement','K');
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


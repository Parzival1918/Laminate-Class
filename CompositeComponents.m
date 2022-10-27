classdef CompositeComponents < handle
    %COMPOSITECOMPONENTS Summary: Use this class to add material properties
    %of the fibre and matrix materials. Then these properties can be used
    %to create a lamina.
    %   Detailed explanation goes here

    properties 
        reinforcementE1 (1,1) {mustBeFloat} = 0 %Reinforcement young's 
                                                % modulus
        matrixE1 (1,1) {mustBeFloat} = 0
        reinforcementv12 (1,1) {mustBeFloat} = 0
        matrixv12 (1,1) {mustBeFloat} = 0
        reinforcementG12 (1,1) {mustBeFloat} = 0
        matrixG12 (1,1) {mustBeFloat} = 0
        reinforcementK (1,1) {mustBeFloat} = 0 %Reinforcement bulk modulus
        matrixK (1,1) {mustBeFloat} = 0

        E_units char {mustBeMember(E_units,{'GPa','MPa','Pa'})} = 'GPa'
        v_units char {mustBeMember(v_units,{'unitless'})} = 'unitless'
        G_units char {mustBeMember(G_units,{'GPa','MPa','Pa'})} = 'GPa'
        K_units char {mustBeMember(K_units,{'GPa','MPa','Pa'})} = 'GPa'
    end
    
    methods
        function obj = CompositeComponents(varargin)
            %COMPOSITECOMPONENTS Construct an instance of this class
            %   Detailed explanation goes here
            narginchk(1,8);
            switch nargin
                case 1
                    obj.reinforcementE1 = varargin{1};
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
    end
end


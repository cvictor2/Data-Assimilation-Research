function obj = override_properties(obj, varargin)
            %Function to override values of existing settings object. It
            %basically runs the constructor again, but is used to override
            %existing settings values.
            %The purpose of this function is to change settings on a
            %restart of a run.

            % !!! WARNING: correct usage is to call the function as
            % follows:
            % obj = override_properties(obj, ...);
            % the following is improper:
            % override_properties(obj, ...);
            % Doing this will not save the changes, except for objects that
            % are passed by reference (the DA_obs 



            % Create the input parser
            p = inputParser;

            % Get a list of all properties of this class
            propertyList = properties(obj);
            
            % For each property, add it as a parameter to the inputParser
            for i = 1:length(propertyList)
                propName = propertyList{i};
                addParameter(p, propName, obj.(propName)); % default is current value
            end

            % Parse the inputs
            parse(p, varargin{:});

            % Update properties based on input
            for i = 1:length(propertyList)
                propName = propertyList{i};
                obj.(propName) = p.Results.(propName);
            end

            
  
        end
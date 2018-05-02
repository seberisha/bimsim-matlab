classdef sphere_bimsim
    % write a description of the class here.
    
    properties(Constant = true)
        DAYS_PER_YEAR =  365;
        MONTHS_PER_YEAR = 12;
        WEEKS_PER_YEAR  = 52;
    end
    
    properties(GetAccess = 'public', SetAccess = 'private')
        % public read access, but private write access.
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        % private read and write access
    end
    properties
        % define the properties of the class here, (like fields of a struct)
        minute = 0;
        hour;
        day;
        month;
        year;
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = mydate(minute,hour,day,month,year)
            % class constructor
            if(nargin > 0)
                obj.minute = minute;
                obj.hour   = hour;
                obj.day    = day;
                obj.month  = month;
                obj.year   = year;
            end
        end
        
        function obj = rollDay(obj,numdays)
            
            obj.day = obj.day + numdays;
        end
        
    end
    
    methods(Access = private)
        function sec = calcSecs(obj)
            sec = obj.minute*60 + obj.hour*60*60 + obj.day*24*60*60;
        end
        
        function TF = isValid(obj)
            TF = obj.minute >= 0 && obj.minute <= 60;
        end
    end
    
    methods(Static = true)
        function printCurrentDate()
            display(datestr(now));
        end
    end
end


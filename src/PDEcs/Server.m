classdef Server < handle
    
    properties (Access = private)
        id_
    end
    
    methods
        
        function this = Server()
            this.id_ = Server_('new');
        end
        
        function delete(this)
            
            Server_('delete', this.id_);
        end
        
        function [ret] = sprod(this, u, v)
            [ret] = Server_('sprod', this.id_, u, v);
        end
        
        function [ret] = pprod(this, u ,v)
            [ret] = Server_('pprod', this.id_, u, v);
        end
    end
    
end


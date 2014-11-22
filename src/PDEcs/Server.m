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
        
        function [ret] = prod(this, u, v)
            [ret] = Server_('prod', this.id_, u, v);
        end
    end
    
end


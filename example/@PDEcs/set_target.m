function set_target(this) 
    % var as true solution
    this.static_var.data  = -(this.static_var.S - this.static_var.true*this.static_var.M)\this.static_var.LoadVector;
end


function f = objective(this, ~)

%             this.set_volatile_var(x);
%             
%             this.volatile_var.u = -(this.static_var.S - this.volatile_var.M)\this.static_var.LoadVector;
    diff = this.volatile_var.u - this.static_var.data;
    % with regularization term
    f = 0.5*diff'* this.static_var.Q*diff + 0.5*this.beta*this.x*this.x;
end


function g = gradient(this, x)

    this.set_volatile_var(x);

    this.volatile_var.u = -(this.static_var.S - this.volatile_var.M)\this.static_var.LoadVector;
    diff = this.volatile_var.u - this.static_var.data;
    this.volatile_var.v = (this.static_var.S - this.volatile_var.M)\(this.static_var.Q*diff);


    g = this.volatile_var.u'* this.static_var.M*this.volatile_var.v...
        + this.beta*this.x;


end

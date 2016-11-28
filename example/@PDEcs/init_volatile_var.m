function init_volatile_var(this)
    % declaration as placehold
    % assuming x0 as constant 1.
    this.volatile_var.M = ...
        this.static_var.fem.assema(this.x0);

    this.x = this.x0;
end

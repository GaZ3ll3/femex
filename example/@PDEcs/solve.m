function solve(this)
    [this.x, this.info] = ipopt(this.x0, this.funcs, this.options);
end

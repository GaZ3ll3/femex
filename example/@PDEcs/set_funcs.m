function set_funcs(this)
    this.funcs.objective = @this.objective;
    this.funcs.gradient  = @this.gradient;
%             this.funcs.iterfunc  = @this.callback;
end


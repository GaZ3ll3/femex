function load_source_functions(this, func_holder) 
    funcs = func_holder();
    for f_id = 1:length(funcs)
        this.sources{f_id} = funcs{f_id};
    end
    
end


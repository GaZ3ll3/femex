function load_source_functions(o, func_holder)
%LOAD_SOURCE_FUNCTIONS loads all sources for qpat.
    funcs = func_holder();
    for f_id = 1:length(funcs)
        o.sources{f_id} = funcs{f_id};
    end
    
end


function monitor(fem)

fields = fieldnames(fem);
len    = size(fields, 1);

sum = 0;

if (len ~= 0) 
    fprintf('-------------------------------------\n');
    fprintf('Name\t\t  Class\t\t Size\n');
    fprintf('-------------------------------------\n\n');
end
for i = 1 : len
    curr = fem.(fields{i});
    fprintf('%s', fields{i});
    
    if isstruct(curr)
        fprintf('\n')
        subfields = fieldnames(curr);
        sublen = size(subfields, 1);
        for j = 1:sublen
            
            str = subfields{j};
            s_len = size(str, 2);
            
            fprintf('    |-%s', str);
            for k = 1: 12 - s_len
                fprintf(' ');
            end
            
            type = class(curr.(str));
            t_len = size(type, 2);
            
            fprintf('%s', type);
            
            for k = 1: 12 - t_len
                fprintf(' ');
            end
            fprintf('\t %d x %d', size(curr.(str), 1), size(curr.(str), 2));
            % get size
            
            if isa(curr.(str), 'double')
                sum = sum + 8 * size(curr.(str), 1) * size(curr.(str), 2);
            elseif isa(curr.(str), 'int32')
                sum = sum + 4 * size(curr.(str), 1) * size(curr.(str), 2);
            end
            
            
            fprintf('\n');
        end
    else
        
        
        s_len = size(fields{i}, 2);
        
        for k = 1: 18 - s_len
            fprintf(' ');
        end
        
        type = class(curr);
        t_len = size(type, 2);
        
        fprintf('%s', type);
            
        for k = 1: 15 - t_len
            fprintf(' ');
        end
        fprintf('%d x %d', size(curr, 1), size(curr, 2));     
        
        
        if isa(curr, 'double')
           sum = sum +  8 * size(curr, 1) * size(curr, 2);
        elseif isa(curr, 'int32')
            sum = sum + 4 * size(curr, 1) * size(curr, 2);
        end
    end
    fprintf('\n');
end

if (len ~= 0) 
    fprintf('\n-------------------------------------\n');
    fprintf('total memory allocated: %d bytes\n', sum);
    fprintf('-------------------------------------\n\n');
end


end
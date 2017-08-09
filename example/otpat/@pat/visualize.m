function visualize(this, var, n, m)
    L = length(var);
    if (L ~= m * n) 
        sprintf('Incorrent dimension of subplot.');
    else 
        figure;
        for i = 1:L
        subplot(n, m, i);
        trisurf(this.model.TriMesh', ...
            this.model.Promoted.nodes(1,:)', ...
            this.model.Promoted.nodes(2,:)', ...
            var{i}, 'EdgeColor', 'none' ); shading interp; colormap jet; colorbar;view(2);
        end
    end
        
    
        
end
        


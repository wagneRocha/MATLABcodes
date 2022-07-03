function treeBranches = vertexOrderBranches(n)
    
    caso = mod(n,5);
    if(caso == 0)
        % nao tem H3 e H2
        desloca = 2;
    elseif(caso == 1)
        % tem H2 ou H3
        desloca = 1;
    elseif(caso == 2)
        % tem H2 e H3
        desloca = 0;
    else
        fprintf('ERRO\n');
    end

    treeBranches = zeros(n,1);
    
    if(caso == 0)
        treeBranches(5) = 1;
    elseif(caso == 1)
        treeBranches(5) = 1;
    else
        treeBranches(6 - desloca) = 1; 
    end
    % DDGP hc Order: {H3, H2, H1, N1, CA1, HA1, C1,
    %                 ...,
    %                 Ni, HZi, CAi, Ci, HAi
    %                 ...,
    %                 Np, HZp, CAp, Cp, HAp}
    % -------------------- 2nd to p-th Residue --------------------
    vertexWithBranches = [8-desloca : 5 : n, 11-desloca : 5 : n];
    treeBranches(vertexWithBranches) = 1;
end
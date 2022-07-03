function tauAngle = vertexOrderTorsionAngles(X, cliques, treeBranches, phibar, psibar)

    n = size(cliques,1);
    tauAngle = nan(n,1);
    
    caso = mod(n,5);
    if(caso == 0)
        % nao tem H3 e H2
        desloca = 0;
    elseif(caso == 1)
        % tem H2 ou H3
        desloca = 1;
    elseif(caso == 2)
        % tem H2 e H3
        desloca = 2;
    else
        fprintf('ERRO\n');
    end
    
    i1 = 1; i2 = 2;
    for i = 4 : n
        if(treeBranches(i))
            if(i == 6)
                
            elseif(mod(i-desloca,5) == 1)
                tauAngle(i) = psibar(i1);
                i1 = i1 + 1;
            elseif(mod(i-desloca,5) == 4)
                tauAngle(i) = phibar(i2);
                i2 = i2 + 1;
            end
        else
            x1 = X(cliques(i,4),:);
            x2 = X(cliques(i,3),:);
            x3 = X(cliques(i,2),:);
            x4 = X(cliques(i,1),:);
            tauAngle(i) = torsionAngleWithPoints(x1, x2, x3, x4)*signTorsionAngle(x1, x2, x3, x4);
        end
    end
end
% I = [atomo_i, [atomo_j (predecessor adjacente), dij, Dij]]
function [I, kv] = instance2MyFormat(I0, m, vertexOrder)

	n = I0(m,1);
	kv = zeros(n,1);
	for i = 1 : n
		kv(i) = length(find(I0(:,1) == i));
	end
	
	I = zeros(n, max(kv), 3);
  
	I(1,:,:) = zeros(max(kv),3);
	
	I(2,:,:) = zeros(max(kv),3);
	I(2,1,:) = [I0(1,2), I0(1,3), I0(1,4)];
	
	I(3,:,:) = zeros(max(kv),3);
	I(3,1,:) = [I0(3,2), I0(3,3), I0(3,4)];
	I(3,2,:) = [I0(2,2), I0(2,3), I0(2,4)];
	
	lastInstanceLineVertexVm1 = 3;
	for i = 4 : n
		orderVertexV = zeros(kv(i),1);
		orderVertexV(1) = vertexOrder(i,2);
		orderVertexV(2) = vertexOrder(i,3);
		orderVertexV(3) = vertexOrder(i,4);
		% ordena os vertices = [DDGP + sobresalentes em orderm decrescente]
		adjacentPredecessorsVertexV = I0(lastInstanceLineVertexVm1 + 1 : lastInstanceLineVertexVm1 + kv(i), 2);
		if(kv(i) > 3)
			DDGPorderVertexV = vertexOrder(i,2:4);
			adjPredVertexV = setdiff(adjacentPredecessorsVertexV, DDGPorderVertexV);
			for l = 1 : kv(i)-3
				orderVertexV(l+3) = adjPredVertexV(kv(i)-3-l+1);
			end
		end
		I0VertexV = zeros(kv(i),3);
		I0VertexVaux = I0(lastInstanceLineVertexVm1 + 1 : lastInstanceLineVertexVm1 + kv(i), 2 : 4);
		for k = 1 : kv(i)
			l = find(adjacentPredecessorsVertexV == orderVertexV(k))';
			I0VertexV(k,:) = I0VertexVaux(l,:);
		end

		j = 1;
		for l = 1 : kv(i)
			I(i,j,:) = [I0VertexV(l,1), I0VertexV(l,2), I0VertexV(l,3)];
			j = j + 1;
		end
		lastInstanceLineVertexVm1 = lastInstanceLineVertexVm1 + kv(i);
	end
end

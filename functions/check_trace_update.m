function res = check_trace_update(A, edges, rob, miobi)
	% Check the computation of the update of the trace of the exponential in dense arithmetic 
	tA = A;
	eA = expm(full(double(A)));
	for j=1:size(edges, 1)
		if strcmp(miobi, 'break')
			tA(edges(j, 1), edges(j, 2)) = 0;
			tA(edges(j, 2), edges(j, 1)) = 0;
		else
			tA(edges(j, 1), edges(j, 2)) = 1;
			tA(edges(j, 2), edges(j, 1)) = 1;
		end
	end
	tA = expm(full(double(tA)));
	res = abs(trace(tA - eA) - rob);
end

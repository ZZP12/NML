function [FLOW,FVAL,UPTAKE,EXITFLAG] = FluxDriver(dataDictionary)

	% Get some stuff from the DF -
	STM = dataDictionary('stoichiometric_matrix');
	OBJVECTOR = dataDictionary('objective_coefficient_array');
	% Get Flux bounds from the DF -
	FluxBounds = dataDictionary('default_flux_bounds_array');
	FluxLB = FluxBounds(:,1);
	FluxUB = FluxBounds(:,2);
	% Get species bounds
	SpeciesBounds = dataDictionary('species_bounds_array');
	% Get numbers
	NUM_Unbalanced = dataDictionary('extra_species_num');
	NUM_Speices = size(SpeciesBounds, 1);
	NUM_Balanced = NUM_Speices - NUM_Unbalanced;

	% Equality constraints
	Aeq = STM((NUM_Unbalanced+1):NUM_Speices, :);
	% Formulate the bV -
	bVEq = zeros(NUM_Balanced,1);

	% Inequality constraints
	UNBALANCED_STM = STM(1:NUM_Unbalanced, :);
	bVLB = SpeciesBounds(1:NUM_Unbalanced, 1);
	bVUB = SpeciesBounds(1:NUM_Unbalanced, 2);
	A = [UNBALANCED_STM ; -1*UNBALANCED_STM];
	bV = [bVUB ; -1*bVLB];

	% Set some values for the options parameter of LINPROG
	options = optimset('TolFun',1e-6);

	% Call the LP solver -
	[FLOW,FVAL,EXITFLAG,OUT,LAM] = linprog(OBJVECTOR, A,bV, Aeq,bVEq, FluxLB,FluxUB, options);

	% Problem2 -- Post-processing - catch the minimal overall flux distribution
    if (EXITFLAG >= 0)
        ObjPro2 = ones([length(OBJVECTOR), 1]);
        A2 = [A; OBJVECTOR];
        bV2 = [bV; FVAL];
        [FLOW,FluxSum,EXITFLAG,OUT2,LAM2] = linprog(ObjPro2, A2,bV2, Aeq,bVEq, FluxLB,FluxUB, options);
    end

	UPTAKE = UNBALANCED_STM*FLOW;
return;

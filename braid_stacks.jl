using Chevie

##################################
# Some functions for convenience #
##################################

subscript(n::Integer) = replace(string(n),
    '-' => '₋',
    '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
    '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉',
)

weyl_wordstring(v::AbstractVector{<:Integer}; sep="") =
    isempty(v) ? "1" : join(("s" * subscript(i) for i in v), sep)

braid_wordstring(v::AbstractVector{<:Integer}; sep="") =
    isempty(v) ? "1" : join(("b" * subscript(i) for i in v), sep)

function is_periodic(G,beta)
	"""
	Check if the braid beta is periodic, i.e., beta^m = full_twist^d, and returns slope d/m.
	
	Input: 
	
	- G: 	a finite reductive group object.
	- beta:	a positive braid in the monoid of positive braids associated to G.
	
	Output:
	
	- A tuple (true, slope) or (false, nothing).
	
	"""
	W = G
	B = beta.M
	full_twist = B(longest(W))^2
	cox_num = div(length(roots(W)), rank(W))

	for m in 1:cox_num
		beta_pow = beta^m
		for d in 1:cox_num
			if gcd(d,m) == 1 && beta_pow == full_twist^d
				return true, d//m
			end
		end
	end
	
    return false, nothing
end

##########################
# Functions for the user #
##########################

function count_points(G,vect,d; double_check=false, output=false, table=true)	
	
	"""
	Compute points on the braid stack M(β,γ) for all unipotent classes γ in G
	
	Input: 
	
	- G: 	a finite reductive group object.
	- vect: a vector of indices of braid group generators which defines the intermediate positive braid. 
	- d: 	an integer exponent applied to the intermediate positive braid defined by vect to obtain the final positive braid.
	
	Options:
	
	- double_check: 	if the braid is periodic, compute Hecke algebra character values directly and check agreement
	- output: 			return the point counts as a vector of pairs (γ,|M(β,γ)^F|)
	- table: 			print human-readable table of the point counts
	
	For example, count_points(G = coxgroup(:G,2), vect = [1,2], d = 3) 
	corresponds to the choice G = G2 and β = (positive lift of s1 * s2)^3.
	"""

	q = big(1) * Mvp(:q)
	W = G
	B = BraidMonoid(W)

	braid = reduce(*,B.(vect))^d
	periodic, slope = is_periodic(G,braid)

	if periodic
		# Compute ϕ_q(braid) using Schur's lemma
		W_ct = CharTable(W)
		W_char_info = charinfo(W)
		W_char_names = W_ct.charnames
		W_char_contents = length(roots(W)) .- (W_char_info.A + W_char_info.a)
		W_char_vals = W_ct.irr[:,position_class(W,image(braid))]
		H_char_vals = (q .^ (slope .* W_char_contents)) .* W_char_vals
	else
		# Compute ϕ_q(braid) directly
		H = hecke(W,q)
		T = Tbasis(H)
		T_braid = T(word(braid))
		H_char_vals = char_values(T_braid)
	end

	if double_check && periodic
		H = hecke(W,q)
		T = Tbasis(H)
		T_braid = T(word(braid))
		true_H_char_vals = char_values(T_braid)
		@assert true_H_char_vals == H_char_vals "H character values incorrectly computed"
		perm = indexin(W_char_names, CharTable(H).charnames)
		@assert true_H_char_vals[perm] == H_char_vals "H character values incorrect even after alignment"
	end
	
	# Count points on M(β,γ) via the formula 
	# |M(β,γ)^F| = (1/|G^F|) * sum_{g in γ^F} sum_{ϕ in Irr(W)} ρ_ϕ(g) * ϕ_q(braid)

	ucl = UnipotentClasses(G)
	uval = UnipotentValues(ucl;classes=true)
	xt = XTable(ucl;classes=true)
	centraliser_sizes_inverted = map(f -> 1//f,xt.centClass)
	rational_unipotent_classes_TeX_names = map(label -> name(TeX(rio();class=label[2]),ucl.classes[label[1]]),xt.classes)
	rational_unipotent_classes_names = fromTeX.(Ref(rio()),rational_unipotent_classes_TeX_names)
	
	stack_counts = Array{Any}(nothing,length(rational_unipotent_classes_names),1)

	for i in 1:length(uval.classes)
		point_count = 0*q
		for j in 1:length(H_char_vals)
			point_count +=  H_char_vals[j] * uval.scalar[j,i]
		end
		point_count *= centraliser_sizes_inverted[i]
		stack_counts[i] = point_count
	end
	
	stack_counts_before_rows_summed = stack_counts
	
	# Determine which rational classes correspond to the same geometric class
    rational_geometric_indices = xt.classes # A list of pairs [n,m] with n counting geometric orbits and m counting rational orbits inside the geometric one  
    class_ids = map(x->x[1],rational_geometric_indices) # Sends [n,m] to n
    duplicated = filter(u -> count(==(u), class_ids) > 1, unique(class_ids)) # Find all labels n that occur more than once (meaning multiple rational classes)
    groups = reverse([ findall(==(u), class_ids) for u in duplicated ]) # Groups of rows which need to be summed
	
	# Now sum the rows appropriately
    for group in groups
        summed_row = sum(stack_counts_before_rows_summed[group, :], dims=1)
        rows_before_summed_row = stack_counts_before_rows_summed[1:group[1]-1,:]
        rows_after_summed_row = stack_counts_before_rows_summed[group[end]+1:end,:]
        stack_counts_before_rows_summed = vcat(rows_before_summed_row, summed_row, rows_after_summed_row)
    end
	
	stack_counts_after_rows_summed = stack_counts_before_rows_summed
	
	if table
		# Find geometric unipotent class names
		rational_label_indices = [ x for g in groups for x in g[2:end] ] 
		keep = setdiff(1:length(rational_unipotent_classes_names), rational_label_indices)
		geometric_unipotent_classes_names = rational_unipotent_classes_names[keep] 
		
		# Make the table
		repr_stack_counts = xrepr.(Ref(rio()),CycPol.(stack_counts_after_rows_summed))
		println("The group is G = ",xrepr(rio(),G))
		println("The braid is β = ",braid_wordstring(word(braid)))
		showtable(repr_stack_counts;
		col_labels=["|M(β,γ)^F|"],
		row_labels=geometric_unipotent_classes_names,
		rows_label="γ ⊆ G"
		)
	end
	
	if output
		return map(i->(ucl.classes[i],stack_counts_after_rows_summed[i]),1:length(ucl.classes))
	end
end

function interval_reps(G,vect,d;table=true)

	"""
	Determines the unipotent classes γ1 and γ2 such that { γ : M(β,γ) non-empty } = [γ1, γ2].
	
	Input: see count_points
	
	Output: a vector of tuples (γ1,γ2)
	
	Options:
	
	- table: print human-readable table of the point counts

	For example, interval_reps(G = coxgroup(:G,2), vect = [1,2], d = 5) returns 
	a vector containing the tuple (A1,G2) because O_\nu = A1 when \nu = 5/6 
	(see Jakob-Yun 2023 Table 2, arXiv:2301.10967v2).
	
	Note: this function should always return a vector containing exactly one tuple.
	
	"""
	
	ucl = UnipotentClasses(G)
	classes = ucl.classes
	P = ucl.orderclasses
	class_count_pairs = count_points(G,vect,d;output=true,table=table)
	non_zero_pairs = filter(pair -> pair[2]!=0,class_count_pairs)
	candidate_interval = map(pair -> pair[1],non_zero_pairs)
	interval_reps = Tuple{UnipotentClass,UnipotentClass}[]
	
	for lower_rep in classes
		for upper_rep in classes
			S = interval(P, ≥, lower_rep, ≤, upper_rep)
			if candidate_interval == S
				push!(interval_reps, (lower_rep,upper_rep))
			end
		end
	end

	return interval_reps
end


function count_points_unique(G,vect,d)

	"""
	Given the class γ1 as in the previous function, return a vector containing the count |M(β,γ1)^F|
	
	"""

	reps = interval_reps(G,vect,d;table=false)

	unique_class = reps[1][1]
	polynomials = Any[]
	
	for pair in count_points(G,vect,d;table=false,output=true)
		if unique_class == pair[1]
			push!(polynomials,pair[2])
		end
	end
	
	return polynomials
end

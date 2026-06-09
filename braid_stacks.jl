using Chevie

# In this file, you will find:
# 1. Some functions for convenience (not for user)
# 2. Functions for the user

####################################################
# 1. Some functions for convenience (not for user) #
####################################################

"""
    subscript(n::Integer)

Return n as a string written with subscript digits and sign.

This is used when displaying Weyl group words and braid words in a compact, human-readable form.

Example: subscript(123) returns the string "₁₂₃".

subscript(-45) returns the string "₋₄₅".

"""

subscript(n::Integer) = replace(string(n),
    '-' => '₋',
    '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
    '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉',
)






"""
    weyl_wordstring(v::AbstractVector{<:Integer}; sep="")

Format a vector of simple-reflection indices as a Weyl group word.

Example: weyl_wordstring([1,2,3]) returns the string "s₁s₂s₃".

weyl_wordstring([1,2,3]; sep="*") returns the string "s₁*s₂*s₃".

"""

weyl_wordstring(v::AbstractVector{<:Integer}; sep="") =
    isempty(v) ? "1" : join(("s" * subscript(i) for i in v), sep)







"""
    braid_wordstring(v::AbstractVector{<:Integer}; sep="")

Format a vector of braid-generator indices as a positive braid word.

The empty vector is displayed as "1". Use sep to insert a separator
between adjacent generators.

Example: braid_wordstring([1,2,3]) returns the string "b₁b₂b₃".

braid_wordstring([1,2,3]; sep="*") returns the string "b₁*b₂*b₃".

"""

braid_wordstring(v::AbstractVector{<:Integer}; sep="") =
    isempty(v) ? "1" : join(("b" * subscript(i) for i in v), sep)







"""
    is_periodic(G, beta)

Check whether beta = β is periodic in the positive braid monoid attached to G.

Note β is an element of the positive braid monoid.
The positive braid monoid is constructed using, for example, the command 
B = BraidMonoid(G)
and a positive braid is constructed using, for example, the command
B([1,2])^3
which is the braid β = b₁b₂b₁b₂b₁b₂.

The function searches for coprime positive integers d and m such that
β^m equals the dth power of the full twist. If successful, it returns
(true, d//m), where d//m is the slope. Otherwise it returns
(false, nothing).

Note we only search up to m = the Coxeter number, because if beta 
is periodic with slope d/m, then m divides the Coxeter number.

Example: for G = G2, the braid β = (positive lift of s1 * s2)^3 
is periodic with slope 1/2, because β^2 = full_twist. This means the function call
is_periodic(coxgroup(:G,2), B(coxgroup(:G,2))([1,2])^3) 
returns (true, 1//2).

On the other hand, the braid β = b₁ is not periodic, so the function call
is_periodic(coxgroup(:G,2), B(coxgroup(:G,2))([1])) 
returns (false, nothing).

"""
function is_periodic(G,beta)
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













#############################
# 2. Functions for the user #
#############################

"""
    count_points(G, vect, d; double_check=false, output=false, table=true)

Compute the point counts |M(β,γ)^F| for each geometric unipotent classes γ in G.

The braid β is defined as the dth power of the positive lift of the word encoded by vect.
Specifically, if vect = [i₁, i₂, ..., iₖ] then β = (b_{i₁} b_{i₂} ... b_{iₖ})^d.

# Arguments
- G: a finite reductive group object.
- vect: indices of braid generators defining the positive braid b_{i₁} b_{i₂} ... b_{iₖ}.
- d: exponent applied to that braid to obtain β = (b_{i₁} b_{i₂} ... b_{iₖ})^d.

# Keywords
- double_check: when true AND the braid is periodic, recompute the Hecke
  algebra character values directly and assert agreement.
- output: when true, return the point counts as (class, count) pairs.
- table: when true, print a human-readable table of the counts.

# Returns
Returns the vector of (class, count) pairs only when output=true;
otherwise the function is used for its printed table.

# Example: The function call

count_points(coxgroup(:G,2), [1,2], 2)

prints the following:

The group is G = G₂
The braid is β = b₁b₂b₁b₂
┌──────┬──────────┐
│γ ⊆ G │|M(β,γ)^F|│
├──────┼──────────┤
│1     │         0│
│A₁    │         0│
│Ã₁    │         0│
│G₂(a₁)│         1│
│G₂    │        x²│
└──────┴──────────┘

"""

function count_points(G,vect,d; double_check=false, output=false, table=true)	
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
	uval = UnipotentValues(ucl;classes=true).scalar
	xt = XTable(ucl;classes=true)
	centraliser_sizes_inverted = map(f -> 1//f,xt.centClass)
	rational_unipotent_classes_TeX_names = map(label -> name(TeX(rio();class=label[2]),ucl.classes[label[1]]),xt.classes)
	rational_unipotent_classes_names = fromTeX.(Ref(rio()),rational_unipotent_classes_TeX_names)
	
	stack_counts = Array{Any}(nothing,length(rational_unipotent_classes_names),1)

	for i in 1:length(xt.classes) # number of rational unipotent classes
		point_count = 0*q
		for j in 1:length(H_char_vals) # number of principal unipotent characters
			point_count +=  H_char_vals[j] * uval[j,i]
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















"""
    interval_reps(G, vect, d; table=true)

Determine interval representatives (γ1, γ2) such that the
unipotent classes γ with non-empty braid stack M(β,γ) form exactly
the interval [γ1, γ2] in the poset of unipotent classes ordered by the closure relation.

The arguments are the same as for count_points.

# Keywords
- table: display the table of point counts for reference (this is helpful when visually checking the result).

# Returns
A tuple (γ1, γ2) of unipotent classes, or nothing if the classes with non-empty braid stack do not form an interval.

Example: The function call

interval_reps(coxgroup(:G,2), [1,2], 2)

returns the tuple (UnipotentClass(G₂(a₁)), UnipotentClass(G₂)) because 
when G = G2 and β = (b₁b₂)^2, the point count |M(β,γ)^F| is non-zero if and only if 
γ is in the interval [G₂(a₁), G₂] in the poset of unipotent classes of G₂.

Another example: The function call

interval_reps(coxgroup(:F,4), [1,2], 5) 

returns the tuple (A1,F4). This example should be contrasted with 
Jakob and Yun's paper arXiv:2301.10967v2. In this paper, 
ν = d/m is a slope which determines a braid β_ν, and
O_ν is alternate notation for the lower representative.
In this case, G = F4 and ν = 5/6 yield the braid β_ν = (b₁b₂)^5.
The function call above returning (A1,F4) agrees with the fact that 
O_ν = A1 for G = F4 and ν = 5/6, as shown in Jakob-Yun's paper 
(see Table 5 of arXiv:2301.10967v2). 

"""

function interval_reps(G,vect,d;table=true)
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

	if length(interval_reps) == 1
		return interval_reps[1]
	else
		println("The conjugacy classes with non-zero point-count do not form an interval")
		return nothing
	end
end


















"""
    count_points_lower(G, vect, d)

Return the point count attached to the lower endpoint selected by interval_reps(G, vect, d).

More precisely, interval_reps(G, vect, d) attempts to find an interval (γ1, γ2) in the poset 
of unipotent classes such that M(β, γ) is non-empty exactly for γ in [γ1, γ2].

This function then extracts the point count |M(β, γ1)^F| attached to the lower representative γ1.

The arguments G, vect, and d determine the braid
β = (b_{i_1} b_{i_2} ... b_{i_r})^d
in the same way as in count_points and interval_reps.

# Returns
The polynomial |M(β,γ1)^F| corresponding to the lower endpoint γ1, provided interval_reps returns 
a unique interval and count_points produces a unique count for that class.

# Errors
Throws an error if interval_reps returns nothing, meaning the classes with
non-zero point count do not form a single interval.

Throws an error if more than one polynomial is found for the lower representative.

# Example: The function call

interval_reps(coxgroup(:G,2), [1,2], 2)

returns the tuple (UnipotentClass(G₂(a₁)), UnipotentClass(G₂)). 

The lower representative is UnipotentClass(G₂(a₁)). 

Then the function call

count_points_lower(coxgroup(:G,2), [1,2], 2)

returns the point-count |M(β, G₂(a₁))^F|, where β = (b₁b₂)^2.

"""

function count_points_lower(G,vect,d)
	reps = interval_reps(G,vect,d;table=false)

	if reps == nothing
		error("No lower representative")
	end

	unique_class = reps[1]
	polynomials = Any[]
	
	for pair in count_points(G,vect,d;table=false,output=true)
		if unique_class == pair[1]
			push!(polynomials,pair[2])
		end
	end
	
	if length(polynomials) == 1
		return polynomials[1]
	else
		error("More than one polynomial/lower representative")
	end
end








"""
    springer_element(G, m)

Returns a reduced word for m-Springer element in the Weyl group of G = coxgroup(:letter,rank)
using the appendix of:  

Broué, Michel; Michel, Jean.  
Sur certains éléments réguliers des groupes de Weyl et les variétés de Deligne-Lusztig associées.(French)  
[Some regular elements of Weyl groups and the associated Deligne-Lusztig varieties]  
Finite reductive groups (Luminy, 1994), 73–139. Progr. Math., 141  





The output of this function is intended to be fed directly into the previous three functions, for example as follows:

count_points(coxgroup(:E,7), springer_element(coxgroup(:E,7), 9), 4)

interval_reps(coxgroup(:E,7), springer_element(coxgroup(:E,7), 9), 4)

count_points_lower(coxgroup(:E,7), springer_element(coxgroup(:E,7), 9), 4)





These replace hand-typing the long and error-prone function calls:

count_points(coxgroup(:E,7), [1,3,4,2,4,7,6,5,4,1,3,4,2,4,7,6,5,4], 4)

interval_reps(coxgroup(:E,7), [1,3,4,2,4,7,6,5,4,1,3,4,2,4,7,6,5,4], 4)

count_points_lower(coxgroup(:E,7), [1,3,4,2,4,7,6,5,4,1,3,4,2,4,7,6,5,4], 4)


"""




function springer_element(G,m)
    m in regular_numbers(G) || error("m must be a regular number for G")
    
    dummy = match(r"^coxgroup\(:([A-Za-z]),(\d+)\)$", string(G))
    dummy === nothing && error("G must look like coxgroup(:A,3)")
    letter, rank = dummy.captures[1], parse(Int, dummy.captures[2])

    if letter == "A"
        floor = fld(rank, 2)
        c = [1:floor; rank:-1:floor+1]
        w = [1:floor; rank:-1:floor]
        if (rank + 1) % m == 0
            i = div(rank + 1, m)
            return vcat(repeat(c, i)...)
        elseif rank % m == 0
            i = div(rank, m)
            return vcat(repeat(w, i)...)
        end
    end


    if letter == "B" || letter == "C"
        c = [1:2:rank; 2:2:rank]
        if (2 * rank) % m == 0
            i = div(2 * rank, m)
            return vcat(repeat(c, i)...)
        end
    end

    if letter == "D"
        c = [3:2:rank; 1; 2:2:rank]
        w = [1:rank; 2:(rank-1)]
        if (2 * rank - 2) % m == 0
            i = div(2 * rank - 2, m)
            return vcat(repeat(c, i)...)
        elseif rank % m == 0
            i = div(rank, m)
            return vcat(repeat(w, i)...)
        end
    end

    if letter == "G"
        c = [1,2]
        if 6 % m == 0
            i = div(6,m)
            return vcat(repeat(c, i)...)
        end
    end

    if letter == "F"
        c = [1,3,2,4]
        w = [1,4,2,3,2,3]
        if 12 % m == 0
            i = div(12, m)
            return vcat(repeat(c, i)...)
        elseif 8 % m == 0
            i = div(8, m)
            return vcat(repeat(w, i)...)  
        end
    end

    if letter == "E" && rank == 6
        c = [1,4,6,2,3,5]
        w = [1,3,4,2,4,6,5,4]
        w_prime = [3,5,4,1,6,3,5,4,2]
        if 12 % m == 0
            i = div(12, m)
            return vcat(repeat(c, i)...)
        elseif 9 % m == 0
            i = div(9, m)
            return vcat(repeat(w, i)...)
        elseif 8 % m == 0
            i = div(8, m)
            return vcat(repeat(w_prime, i)...)
        end
    end

    if letter == "E" && rank == 7
        c = [1,4,6,2,3,5,7]
        w = [1,3,4,2,4,7,6,5,4]
        if 18 % m == 0
            i = div(18, m)
            return vcat(repeat(c, i)...)
        elseif 14 % m == 0
            i = div(14, m)
            return vcat(repeat(w, i)...)
        end 
    end

    if letter == "E" && rank == 8
        c = [1,4,6,8,2,3,5,7]
        w = [1,3,4,2,4,8,7,6,5,4]
        w_prime = [8,7,6,5,4,2,3,1,4,3,5,4]
        if 30 % m == 0
            i = div(30,m)
            return vcat(repeat(c,i)...)
        elseif 24 % m == 0
            i = div(24,m)
            return vcat(repeat(w,i)...)
        elseif 20 % m == 0
            i = div(20,m)
            return vcat(repeat(w_prime,i)...)
        end
    end
end
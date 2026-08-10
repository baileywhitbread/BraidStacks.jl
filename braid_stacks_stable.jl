using Chevie

# In this file, you will find:

# 1. Functions not for the user
#	1.1 subscript
#	1.2 weyl_wordstring
#	1.3 braid_wordstring
#	1.4 is_periodic

# 2. Functions for the user
#	2.1 count_points
#	2.2 interval_reps
#	2.3 count_points_lower
#	2.4 springer_element
#	2.5 regular_numbers
#	2.6 regular_elliptic_numbers

# 1. Functions not for the user.

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

The empty vector is displayed as "1" (the identity). 

Example: weyl_wordstring([1,2,3]) returns the string "s₁s₂s₃".

Use sep to insert a separator between adjacent generators:

weyl_wordstring([1,2,3]; sep="*") returns the string "s₁*s₂*s₃".

"""

weyl_wordstring(v::AbstractVector{<:Integer}; sep="") =
    isempty(v) ? "1" : join(("s" * subscript(i) for i in v), sep)

"""
    braid_wordstring(v::AbstractVector{<:Integer}; sep="")

Format a vector of braid-generator indices as a positive braid word.

The empty vector is displayed as "1" (the identity). 

Example: braid_wordstring([1,2,3]) returns the string "b₁b₂b₃".

Use sep to insert a separator between adjacent generators:

braid_wordstring([1,2,3]; sep="*") returns the string "b₁*b₂*b₃".

"""

braid_wordstring(v::AbstractVector{<:Integer}; sep="") =
    isempty(v) ? "1" : join(("b" * subscript(i) for i in v), sep)

"""
    is_periodic(G, beta)

Determine whether the positive braid `beta` is periodic in the braid monoid
attached to `G`.

A positive braid is periodic when there are relatively prime integers
`d >= 0` and `m >= 1` such that

    beta^m == full_twist^d.

Since the braid relations preserve word length, the only possible slope is

    d//m = length(word(beta)) // length(word(full_twist)).

The function computes this slope and verifies the corresponding braid
equality. It returns `(true, d//m)` when the equality holds and
`(false, nothing)` otherwise. In particular, the identity braid is periodic
with slope `0//1`.

In rank one, every positive braid is a power of the unique braid generator
and is therefore periodic. This case is handled directly because Chevie can
store equivalent `A₁` braids in different internal forms.

# Examples

```julia
julia> G = coxgroup(:G, 2);

julia> B = BraidMonoid(G);

julia> is_periodic(G, B(1, 2)^3)
(true, 1//2)

julia> is_periodic(G, B(1))
(false, nothing)

julia> is_periodic(G, B())
(true, 0//1)
```

"""

function is_periodic(G, beta)
    B = beta.M
    full_twist = B(longest(G))^2

    beta_length = length(word(beta))
    full_twist_length = length(word(full_twist))

    slope = beta_length // full_twist_length
    d = numerator(slope)
    m = denominator(slope)

    if semisimplerank(G) == 1
        return true, slope
    end

    if beta^m == full_twist^d
        return true, slope
    end

    return false, nothing
end

# 2. Functions for the user.

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
  algebra character values directly (ie. without Trinh 2021 'From the Hecke 
  category to the unipotent locus' Corollary 9.2.2.) and assert agreement.
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

	braid = B(vect...)^d
	periodic, slope = is_periodic(G,braid)
	W_ct = CharTable(W)
	W_char_names = W_ct.charnames

	if periodic
		# Compute ϕ_q(braid) using Trinh 2021 Corollary 9.2.2.
		W_char_info = charinfo(W)
		W_char_contents = length(roots(W)) .- (W_char_info.A + W_char_info.a)
		W_char_vals = W_ct.irr[:,position_class(W,image(braid))]
		H_char_vals = (q .^ (slope .* W_char_contents)) .* W_char_vals
	else
		# Compute ϕ_q(braid) directly
		H = hecke(W,q)
		direct_H_char_vals = char_values(H,word(braid))
		positions = indexin(W_char_names,CharTable(H).charnames)
		any(isnothing,positions) && error("Could not align the Hecke and Weyl character tables")
		H_char_vals = direct_H_char_vals[something.(positions)]
	end

	if double_check && periodic
		# Compute ϕ_q(braid) directly and check with Trinh 2021 Corollary 9.2.2.
		H = hecke(W,q)
		true_H_char_vals = char_values(H,word(braid))
		positions = indexin(W_char_names,CharTable(H).charnames)
		any(isnothing,positions) && error("Could not align the Hecke and Weyl character tables")
		@assert true_H_char_vals[something.(positions)] == H_char_vals "H character values incorrectly computed"
	end
	
	# Count points on M(β,γ) via the formula 
	# |M(β,γ)^F| = (1/|G^F|) * sum_{g in γ^F} sum_{ϕ in Irr(W)} ρ_ϕ(g) * ϕ_q(braid)
	
	# We need the unipotent character table for G^F 
	# UnipotentValues forces us to use Int64 which overflows for groups like B9
	# So we have a small work-around
	
	# Work-around begins
	ucl = UnipotentClasses(G)
	values = UnipotentValues(ucl;q=q,classes=true)
	uval = values.scalar
	# Work-around ends
	
	principal_rows = charnumbers(UnipotentCharacters(W).harishChandra[1])
	length(principal_rows) == length(H_char_vals) ||
		throw(DimensionMismatch("principal-series and Hecke character counts differ"))
	centraliser_sizes_inverted = map(f -> 1//f,values.centClass)
	
	geometric_counts = Any[zero(q) for _ in ucl.classes]

	for i in 1:length(values.classes) # number of rational unipotent classes
		point_count = 0*q
		for (j,row) in enumerate(principal_rows) # principal unipotent characters
			point_count +=  H_char_vals[j] * uval[row,i]
		end
		point_count *= centraliser_sizes_inverted[i]
		geometric_index = first(values.classes[i])
		geometric_counts[geometric_index] += point_count
	end
	
	geometric_counts = improve_type.(geometric_counts)
	
	if table
		# Find geometric unipotent class names
		geometric_unipotent_classes_names = name.(Ref(rio()),ucl.classes)
		
		# Make the table
		repr_stack_counts = reshape(xrepr.(Ref(rio()),CycPol.(geometric_counts)),:,1)
		println("The group is G = ",xrepr(rio(),G))
		println("The braid is β = ",braid_wordstring(word(braid)))
		showtable(repr_stack_counts;
		col_labels=["|M(β,γ)^F|"],
		row_labels=geometric_unipotent_classes_names,
		rows_label="γ ⊆ G"
		)
	end
	
	if output
		return map(i->(ucl.classes[i],geometric_counts[i]),1:length(ucl.classes))
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

The function is written in two parts (interval_reps and _interval_reps) 
to speed up the next function, count_points_lower

"""

function interval_reps(G,vect,d;table=true)
	ucl = UnipotentClasses(G)
	P = ucl.orderclasses
	class_count_pairs = count_points(G,vect,d;output=true,table=table)
	reps = _interval_reps(P,class_count_pairs)

	return reps
end

function _interval_reps(P,class_count_pairs)
	candidate_interval = [class for (class,count) in class_count_pairs if !iszero(count)]
	lower_candidates = minima(P,candidate_interval)
	upper_candidates = maxima(P,candidate_interval)

	if length(lower_candidates) == 1 && length(upper_candidates) == 1
		lower_rep = only(lower_candidates)
		upper_rep = only(upper_candidates)
		S = interval(P, ≥, lower_rep, ≤, upper_rep)
		if issetequal(candidate_interval,S)
			return (lower_rep,upper_rep)
		end
	end

	return nothing
end

"""
    count_points_lower(G, vect, d)

Return the point count attached to the lower endpoint returned by interval_reps(G, vect, d).

More precisely, interval_reps(G, vect, d) attempts to find an interval (γ1, γ2) in the poset 
of unipotent classes such that M(β, γ) is non-empty exactly for γ in [γ1, γ2].

This function then extracts the point count |M(β, γ1)^F| attached to the lower representative γ1.

The arguments are the same as for count_points and interval_reps.

# Returns
The polynomial |M(β,γ1)^F| corresponding to the lower endpoint γ1, provided interval_reps returns 
a unique interval and count_points produces a unique count for that class.

# Errors
Throws an error if interval_reps returns nothing, meaning the classes with
non-zero point count do not form a single interval. (This is actually possible!)

Throws an error unless exactly one point count is found for the lower representative. (Never observed.)

# Example: The function call

interval_reps(coxgroup(:G,2), [1,2], 2)

returns the tuple (UnipotentClass(G₂(a₁)), UnipotentClass(G₂)). 

The lower representative is UnipotentClass(G₂(a₁)). 

Then the function call

count_points_lower(coxgroup(:G,2), [1,2], 2)

returns the point-count |M(β, G₂(a₁))^F|, where β = (b₁b₂)^2.

"""

function count_points_lower(G,vect,d)
	ucl = UnipotentClasses(G)
	P = ucl.orderclasses
	class_count_pairs = count_points(G,vect,d;table=false,output=true)
	reps = _interval_reps(P,class_count_pairs)

	if isnothing(reps)
		println("The unipotent classes with non-zero point counts do not form an interval")
		error("No lower representative")
	end

	lower_class, _ = reps
	lower_counts = (count for (class,count) in class_count_pairs if class == lower_class)
	return only(lower_counts)
end

"""
    springer_element(G, m)

Return a reduced word for an `m`-Springer element in the Weyl group `G`.

The group `G` must be an irreducible crystallographic Weyl group of type
A, B/C, D, E, F, or G. The parameter `m` must be an integer regular number
satisfying `m >= 2`, following Springer's convention.

The returned vector contains simple-reflection indices. Its positive braid
lift is an `m`-th root of the full twist. The constructions follow Appendix 1 of:

Broué, Michel; Michel, Jean.
*Sur certains éléments réguliers des groupes de Weyl et les variétés de
Deligne-Lusztig associées.*

# Example

```julia
julia> springer_element(coxgroup(:G, 2), 3)
4-element Vector{Int64}:
 1
 2
 1
 2
```
"""

function springer_element(G, m)
    components = refltype(G)
    length(components) == 1 ||
        throw(ArgumentError("G must be an irreducible Weyl group"))

    weyl_type = only(components)
    series = weyl_type.series
    supported_series = (:A, :B, :D, :E, :F, :G)
    series in supported_series ||
        throw(ArgumentError("Springer elements are not implemented for type $series"))

    rank = length(indices(weyl_type))

    m isa Integer ||
        throw(ArgumentError("m must be an integer"))
    m >= 2 ||
        throw(ArgumentError("m must satisfy m >= 2"))
    m in regular_numbers(G) ||
        throw(ArgumentError("m = $m is not a regular number for G"))

    candidates = if series == :A
        half_rank = fld(rank, 2)
        coxeter_word = [1:half_rank; rank:-1:half_rank+1]
        period_rank_word = [1:half_rank; rank:-1:half_rank]
        ((rank + 1, coxeter_word), (rank, period_rank_word))
    elseif series == :B
        coxeter_word = [1:2:rank; 2:2:rank]
        ((2 * rank, coxeter_word),)
    elseif series == :D
        coxeter_word = [3:2:rank; 1; 2:2:rank]
        period_rank_word = [1:rank; 2:rank-1]
        ((2 * rank - 2, coxeter_word), (rank, period_rank_word))
    elseif series == :G
        ((6, [1, 2]),)
    elseif series == :F
        ((12, [1, 3, 2, 4]),
         (8, [1, 4, 2, 3, 2, 3]))
    elseif rank == 6
        ((12, [1, 4, 6, 2, 3, 5]),
         (9, [1, 3, 4, 2, 4, 6, 5, 4]),
         (8, [3, 5, 4, 1, 6, 3, 5, 4, 2]))
    elseif rank == 7
        ((18, [1, 4, 6, 2, 3, 5, 7]),
         (14, [1, 3, 4, 2, 4, 7, 6, 5, 4]))
    elseif rank == 8
        ((30, [1, 4, 6, 8, 2, 3, 5, 7]),
         (24, [1, 3, 4, 2, 4, 8, 7, 6, 5, 4]),
         (20, [8, 7, 6, 5, 4, 2, 3, 1, 4, 3, 5, 4]))
    else
        throw(ArgumentError("Springer elements are not implemented for type E$rank"))
    end

    for (period, word) in candidates
        if period % m == 0
            return repeat(word, period ÷ m)
        end
    end

    error("No Springer element construction is available for type $(series)$(rank) and m = $m")
end

"""
    regular_numbers(G)

Return the nontrivial regular numbers for `G`, following Springer's
convention that a regular number satisfies `d >= 2`.

Chevie's `regular_eigenvalues(G)` returns the regular roots of unity.
This function returns their distinct multiplicative orders, excluding
order `1`, in increasing order.

"""

function regular_numbers(G)
    orders = unique(order.(regular_eigenvalues(G)))
    return sort(filter(d -> d >= 2, orders))
end

"""
    regular_elliptic_numbers(G)

Return the nontrivial regular elliptic numbers for `G`, in increasing
order.

A regular number `d` is elliptic when a corresponding regular element
has no nonzero fixed vectors in the reflection representation.
Equivalently, `d` divides none of the exponents `degree - 1` of `G`.
"""
function regular_elliptic_numbers(G)
    exponents = degrees(G) .- 1

    return filter(regular_numbers(G)) do d
        all(e -> e % d != 0, exponents)
    end
end

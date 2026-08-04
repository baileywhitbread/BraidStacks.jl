using Chevie

# This script proves the claims made in
# 'The isoclinic Deligne--Simpson problem and rigid connections'
# by Masoud Kamgarpour and Bailey Whitbread

# One should begin by first loading the functions found at
# https://github.com/baileywhitbread/BraidStacks.jl/raw/main/braid_stacks_stable.jl

# In this file, we prove:
# 1. The correctness of Table 1 in Theorem 1.
# 2. The physical rigidity of the connections in Theorem 2 (a).

# 1. The correctness of Table 1 in Theorem 1.
# Table 1 lists C_ν when G is exceptional and ν = d/m is a slope in lowest terms between 0 and 1 with regular denominator.
# C_ν is the lower class returned by interval_reps, calculated below

for G in [coxgroup(:G,2),coxgroup(:F,4),coxgroup(:E,6),coxgroup(:E,7),coxgroup(:E,8),]
	for m in regular_numbers(G)
		for d in 1:(m-1)
			if gcd(d,m) == 1
				xprint("G = ",G,", ν = ",d//m,": C_ν = ")
				xprintln(first(interval_reps(G,springer_element(G,m),d;table=false)))
			end
		end
	end
end

# 2. The physical rigidity of the connections in Theorem 2 (a).
# Theorem 2 (b) states that if (G,ν) is a pair from Theorem 2 (a) and rank(G) <= 8 then 
# the isoclinic connection with local data (ν,C_ν) is physically rigid, i.e., |M(β_ν,C_ν)^F| = 1. 

groups_rank_eight_or_less = FiniteCoxeterGroup[
    [coxgroup(:A, r) for r in 1:8]...,
    [coxgroup(:B, r) for r in 2:8]...,
    [coxgroup(:C, r) for r in 3:8]...,
    [coxgroup(:D, r) for r in 4:8]...,
    coxgroup(:G, 2),
    coxgroup(:F, 4),
    [coxgroup(:E, r) for r in 6:8]...,
]

for G in groups_rank_eight_or_less
	reg_ell_nums = regular_elliptic_numbers(G)
	
	## Case (i)
	for m in reg_ell_nums
		xprint("G = ",G,", ν = ",1//m,": |M(β_ν,C_ν)^F| = ")
		xprintln(count_points_lower(G,springer_element(G,m),1))
	end
	
	## Case (ii)
	h = maximum(regular_numbers(G))
	xprint("G = ",G,", ν = ",(h+1)//h,": |M(β_ν,C_ν)^F| = ")
	xprintln(count_points_lower(G,springer_element(G,h),h+1))
	
	## Case (iii)
	type = only(refltype(G))
    series = type.series
    rank = type.rank
    slopes = Tuple{Int,Int}[]

    if series == :A
        n = rank + 1
        m = n

        slopes = [
            (d, m) for d in 2:(m - 1)
            if gcd(d, m) == 1 &&
               ((n - 1) % d == 0 || (n + 1) % d == 0)
        ]

    elseif series == :B && type.cartanType == 2
        n = rank
        m = 2 * n

        slopes = [
            (d, m) for d in 2:(m - 1)
            if gcd(d, m) == 1 &&
               ((n + 1) % d == 0 || (2 * n + 1) % d == 0)
        ]

        if iseven(n) && 3 < n && gcd(3, n) == 1
            push!(slopes, (3, n))
        end

    elseif series == :B && type.cartanType == 1
        n = rank
        m = 2 * n

        slopes = [
            (d, m) for d in 2:(m - 1)
            if gcd(d, m) == 1 &&
               ((m - 1) % d == 0 || (m + 1) % d == 0)
        ]

    elseif series == :D
        n = rank
        m = 2 * n - 2

        slopes = [
            (d, m) for d in 2:(m - 1)
            if gcd(d, m) == 1 &&
               ((2 * n) % d == 0 || (2 * n - 1) % d == 0)
        ]

        if iseven(n) && 3 < n && gcd(3, n) == 1
            push!(slopes, (3, n))
        end

    elseif series == :G
        slopes = [(2, 3)]

    elseif series == :F
        slopes = [(3, 8), (5, 8), (3, 4)]

    elseif series == :E && rank == 6
        slopes = [(2, 9), (4, 9), (7, 9)]

    elseif series == :E && rank == 7
        slopes = [(3, 14), (9, 14), (11, 14), (7, 18)]

    elseif series == :E && rank == 8
        slopes = [(3, 10), (2, 15), (4, 15), (7, 15), (3, 20), (13, 20), (19, 24)]
    end

    for (d, m) in slopes
        xprint("G = ",G,", ν = ",d//m,": |M(β_ν,C_ν)^F| = ")
        xprintln(count_points_lower(G,springer_element(G,m),d))
    end
end

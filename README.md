# BraidStacks.jl

Tools to compute the number of points on braid stacks. We heavily rely on Jean Michel's port of the computer algebra system [Chevie](https://github.com/jmichel7/Chevie.jl). 

Tools offered:

- `count_points(G,vect,d)` 

Compute points on the braid stack $`M(β,γ)`$ for all unipotent classes $`γ`$ in $`G`$. Here, $`β`$ is defined by `(G,vect,d)` in the following manner.

Input: 
	
	- `G`: 	a finite reductive group object.
	- `vect`: a vector of indices of braid group generators which defines the intermediate positive braid. 
	- `d`: 	an integer exponent applied to the intermediate positive braid defined by vect to obtain the final positive braid.
	
Options:
	
	- `double_check`: if the braid is periodic, compute Hecke algebra character values directly and check agreement
	- `output`: return the point counts as a vector of pairs (γ,|M(β,γ)^F|)
	- `table`: print human-readable table of the point counts
	
For example, `count_points(G = coxgroup(:G,2), vect = [1,2], d = 3)` corresponds to the choice $`G=G_2`$ and $`β= (\widetilde{s_1s_2})^3`$.

- `interval_reps(G,vect,d)` 

Determines the unipotent classes $`γ_1`$ and $`γ_2`$ such that $`\{ γ \ \text{such that}\ M(β,γ)\ \text{is non-empty} \} = [γ_1, γ_2]`$ in the poset of unipotent orbits.

- `count_points_unique(G,vect,d)` 

Given the class $`γ_1`$ as in the previous function, return a vector containing the count $`|M(β,γ_1)^F|`$.


## Getting started

Download and install [Julia](https://julialang.org/downloads/), and then install [Chevie](https://github.com/jmichel7/Chevie.jl). Next, place the braid_stacks.jl file in Julia's bin folder. In the REPL (Julia's interactive command-line), copy-paste and run the below:

```julia
include("braid_stacks.jl")
```



## Examples

```julia
julia> include("braid_stacks.jl")
count_points_unique (generic function with 1 method)

julia> G=coxgroup(:G,2)
G₂

julia> count_points(G,[1,2],3)
The group is G = G₂
The braid is β = b₁b₂b₁b₂b₁b₂
┌──────┬──────────┐
│γ ⊆ G │|M(β,γ)^F|│
├──────┼──────────┤
│1     │         0│
│A₁    │         0│
│Ã₁    │         1│
│G₂(a₁)│        x²│
│G₂    │        x⁴│
└──────┴──────────┘

julia> interval_reps(G,[1,2],3)
The group is G = G₂
The braid is β = b₁b₂b₁b₂b₁b₂
┌──────┬──────────┐
│γ ⊆ G │|M(β,γ)^F|│
├──────┼──────────┤
│1     │         0│
│A₁    │         0│
│Ã₁    │         1│
│G₂(a₁)│        x²│
│G₂    │        x⁴│
└──────┴──────────┘
1-element Vector{Tuple{UnipotentClass, UnipotentClass}}:
 (UnipotentClass(Ã₁), UnipotentClass(G₂))

julia> count_points_unique(G,[1,2],3)
1-element Vector{Any}:
 1
```

Below is an example for the slope $`\nu=d/m=5/14`$ and the group $`G=E_7`$. In this case, the braid is $`β= (\widetilde{s_4s_2\mathfrak{c}^{-1}})^5`$.

```julia
julia> G=coxgroup(:E,7)
E₇

julia> count_points(G,[4,2,7,6,5,4,3,2,1],5)
The group is G = E₇
The braid is β = b₄b₂b₃b₅b₄b₃b₁b₆b₅b₄b₂b₃b₁b₄b₃b₅b₄b₂b₆b₅b₇b₆b₅b₄b₂b₃b₁b₄b₃b₅b₄b₂b₆b₅b₄b₃b₁b₇b₆b₅b₄b₂b₃b₄b₅
┌─────────┬──────────────────┐
│γ ⊆ G    │        |M(β,γ)^F|│
├─────────┼──────────────────┤
│E₇       │               x³⁸│
│E₇(a₁)   │               x³⁶│
│E₇(a₂)   │       x²⁶Φ₁Φ₂Φ₄Φ₈│
│E₆       │           x²⁶Φ₄Φ₈│
│E₇(a₃)   │       x²⁴Φ₁Φ₂Φ₄Φ₈│
│E₆(a₁)   │    (x⁸+x⁴-1)x²⁰Φ₄│
│D₆       │           x²⁴Φ₄Φ₈│
│E₇(a₄)   │       x²⁰Φ₁Φ₂Φ₄Φ₈│
│A₆       │           x²⁰Φ₄Φ₈│
│D₅+A₁    │    (x⁶+x²-1)x¹⁸Φ₄│
│D₆(a₁)   │    (x⁶+x²-1)x¹⁸Φ₄│
│E₇(a₅)   │(x⁶+x²-1)x¹⁴Φ₁Φ₂Φ₄│
│D₅       │    (x⁴+x²+2)x¹⁸Φ₄│
│D₆(a₂)   │(x⁴+x²+2)x¹⁴Φ₁Φ₂Φ₄│
│E₆(a₃)   │(x⁴+x²+2)x¹⁴Φ₁Φ₂Φ₄│
│A₅′      │    (x⁴+x²+2)x¹⁴Φ₄│
│D₅(a₁)+A₁│(x⁸+x⁴-2x²+1)x¹⁰Φ₄│
│A₅+A₁    │    (x⁶+x²-1)x¹²Φ₄│
│A₄+A₂    │    (x⁶+x²-1)x¹⁰Φ₄│
│A₅″      │             x¹²Φ₄│
│D₅(a₁)   │ (x⁶+x⁴+x²-2)x¹⁰Φ₄│
│D₄+A₁    │    (x⁶+x⁴+x²-1)x⁸│
│A₄+A₁    │(x⁶+2x⁴+2x²-1)x⁸Φ₄│
│D₄       │                x⁸│
│A₃+A₂+A₁ │      x⁴Φ₁Φ₂Φ₃Φ₄Φ₆│
│A₄       │             x⁸Φ₄²│
│A₃+A₂    │             x⁶Φ₄²│
│D₄(a₁)+A₁│                x⁸│
│D₄(a₁)   │          x²Φ₁Φ₂Φ₄│
│A₃+2A₁   │              x⁴Φ₄│
│(A₃+A₁)′ │              x²Φ₄│
│(A₃+A₁)″ │                 0│
│2A₂+A₁   │                Φ₄│
│2A₂      │                 0│
│A₃       │                 0│
│A₂+3A₁   │                 0│
│A₂+2A₁   │                 0│
│A₂+A₁    │                 0│
│4A₁      │                 0│
│A₂       │                 0│
│3A₁″     │                 0│
│3A₁′     │                 0│
│2A₁      │                 0│
│A₁       │                 0│
│1        │                 0│
└─────────┴──────────────────┘

julia> interval_reps(G,[4,2,7,6,5,4,3,2,1],5;table=false)
1-element Vector{Tuple{UnipotentClass, UnipotentClass}}:
 (UnipotentClass(2A₂+A₁), UnipotentClass(E₇))

julia> count_points_unique(G,[4,2,7,6,5,4,3,2,1],5)
1-element Vector{Any}:
 q²+1

julia>
```


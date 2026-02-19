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
	
	For example, `count_points(G = coxgroup(:G,2), vect = [1,2], d = 3)`
	corresponds to the choice $`G=G_2`$ and $`β= (\widetilde{s_1s_2})^3`$.

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

julia>
```


# BraidStacks.jl

Tools to compute the number of points on braid stacks. We heavily rely on Jean Michel's port of the computer algebra system [Chevie](https://github.com/jmichel7/Chevie.jl). 

Tools offered:

- `count_points` 

Compute points on the braid stack $`M(β,γ)`$ for all unipotent classes $`γ`$ in $`G`$.

- `interval_reps` 

Determines the unipotent classes $`γ1`$ and $`γ2`$ such that $`{ γ : M(β,γ) non-empty } = [γ1, γ2]`$.

- `count_points_unique` 

Given the class $`γ1`$ as in the previous function, return a vector containing the count $`|M(β,γ1)^F|`$.


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


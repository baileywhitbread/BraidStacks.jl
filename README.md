# BraidStacks.jl

Tools to compute the number of points on braid stacks and address the isoclinic Deligne--Simpson problem. 
This was written for the paper [arXiv:2603.20499](https://arxiv.org/abs/2603.20499). 
We heavily rely on Jean Michel's port of the computer algebra system [Chevie.jl](https://github.com/jmichel7/Chevie.jl).

## Getting started

1. Download and install [Julia](https://julialang.org/downloads/). 
2. Install [Chevie](https://github.com/jmichel7/Chevie.jl):

```julia
using Pkg; Pkg.add("Chevie")
```

Make sure you are using Chevie v0.1.14 or later.

Check:

```julia
using Pkg; Pkg.status("Chevie")
```

Update Chevie:

```julia
using Pkg; Pkg.update("Chevie")
```

3. Paste the contents of `braid_stacks.jl` into Julia's command line.  
Alternatively, place the file `braid_stacks.jl` in Julia's working directory, then run:

```julia
include("braid_stacks.jl")
```

## Main functions

- `count_points(G, vect, d; output=false, table=true)` computes the point count for every geometric unipotent class. It prints the table when `table=true` and returns the class-count pairs when `output=true`.
- `interval_reps(G, vect, d; table=true)` checks whether the classes with nonzero point count form a single interval and, if so, returns its lower and upper representatives; otherwise it returns `nothing`.
- `count_points_lower(G, vect, d)` returns the point count at that lower representative and throws an error if the nonzero support is not a single interval.
- `springer_element(G, m)` returns a reduced word for an `m`-Springer element. Here `G` must be an irreducible crystallographic Weyl group and `m >= 2` must belong to `regular_numbers(G)`.
- `regular_numbers(G)` returns the nontrivial regular numbers of `G`.
- `regular_elliptic_numbers(G)` returns the nontrivial regular elliptic numbers of `G`.

## Warm up example

Goal: Count points on $$M(β,γ)$$ for all geometric unipotent classes $\gamma$ when $G = G_2$ and $β = (b_1b_2)^2$.

We make the function call

```julia
count_points(coxgroup(:G,2),[1,2],2)
```

This prints the following. Chevie displays the field-size variable `q` as `x` in these tables, and `Φₙ` denotes the $n$th cyclotomic polynomial in `x`.

```julia
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
```


## Using the code for the isoclinic Deligne--Simpson problem

We solve the isoclinic Deligne--Simpson problem by counting points on the braid stack $M(β,γ)$. 

Each slope $\nu$ gives rise to a braid $\beta_\nu$. 

Write $\nu = d/m$ in lowest terms with $m$ a regular number for $W$.

Up to cyclic shift, the braid looks like $$\beta_\nu = \widetilde{w}^d$$.

Here, $w\in W$ is an $m$-Springer element and $\widetilde{w}$ is the positive lift to the braid monoid.

Goal: Count points on $$M(\beta_\nu,\gamma)$$ for all geometric unipotent classes $\gamma$ when $G = F_4$ and $\nu = 5/8$.

In the Springer-element tables below, you will find that the word $w = 142323$ is an $m$-Springer element when $G = F_4$ and $m = 8$.

This means we count points when $\beta_\nu = (142323)^5$. So we call the function

```julia
count_points(coxgroup(:F,4),[1,4,2,3,2,3],5)
```

This prints the following:

```julia
The group is G = F₄
The braid is β = b₁b₂b₁b₃b₂b₁b₃b₂b₃b₄b₃b₂b₁b₃b₂b₃b₄b₃b₂b₁b₃b₂b₃b₄b₁b₂b₄b₃b₂b₃
┌──────┬───────────────┐
│γ ⊆ G │     |M(β,γ)^F|│
├──────┼───────────────┤
│1     │              0│
│A₁    │              0│
│Ã₁    │              1│
│A₁+Ã₁ │          x²Φ₄²│
│Ã₂    │           x⁶Φ₄│
│A₂    │           x⁶Φ₄│
│A₂+Ã₁ │  (x⁶+x²-1)x⁴Φ₄│
│Ã₂+A₁ │  (x⁶+x²-1)x⁶Φ₄│
│B₂    │(x⁶+x⁴+2x²+1)x⁸│
│C₃(a₁)│  (x⁶+x²-1)x⁸Φ₄│
│F₄(a₃)│    x¹⁰Φ₁Φ₂Φ₄Φ₈│
│C₃    │        x¹⁴Φ₄Φ₈│
│B₃    │        x¹⁴Φ₄Φ₈│
│F₄(a₂)│    x¹⁴Φ₁Φ₂Φ₄Φ₈│
│F₄(a₁)│            x²⁴│
│F₄    │            x²⁶│
└──────┴───────────────┘
```

One visually inspects the table and finds $M(\beta_\nu,\gamma)\neq \emptyset$ if and only if $\gamma\geq \widetilde{A_1}$.   

(This is often difficult to determine, unless you have the poset memorised!)  

To check whether the classes with nonzero point count form a single interval and, when they do, obtain its lower and upper representatives, use

```julia
interval_reps(coxgroup(:F,4),[1,4,2,3,2,3],5)
```

This prints the table above and returns the pair whose lower and upper representatives are $\widetilde{A_1}$ and $F_4$, respectively. The precise display of Chevie objects can depend on the output mode.

To count points at this lower representative $C_\nu = \widetilde{A_1}$, use

```julia
count_points_lower(coxgroup(:F,4),[1,4,2,3,2,3],5)
```

which returns the point-count in the row $\gamma = C_\nu = \widetilde{A_1}$.


## Tables of Springer elements

For reference, we give tables of $m$-Springer elements $w$ for each nontrivial regular number $m \geq 2$ of each exceptional group $G$, using the appendix of:

>Broué, Michel; Michel, Jean.  
Sur certains éléments réguliers des groupes de Weyl et les variétés de Deligne-Lusztig associées.(French)  
[Some regular elements of Weyl groups and the associated Deligne-Lusztig varieties]  
Finite reductive groups (Luminy, 1994), 73–139. Progr. Math., 141  

```
G = G_2
m | w
6 | 12 = c
3 | c^2
2 | c^3
```

```
G = F_4
m  | w
12 | 1324 = c
8  | 142323
6  | c^2
4  | c^3
3  | c^4
2  | c^6
```

```
G = E_6
m  | w
12 | 146235 = c
9  | 13424654
8  | 354163542
6  | c^2
4  | c^3
3  | c^4
2  | c^6
```

```
G = E_7
m  | w
18 | 1462357 = c
14 | 134247654 = x
9  | c^2
7  | x^2
6  | c^3
3  | c^6
2  | c^9
```

```
G = E_8
m  | w
30 | 14682357 = c
24 | 1342487654 = x
20 | 876542314354
15 | c^2
12 | x^2
10 | c^3
8  | x^3
6  | c^5
5  | c^6
4  | x^6
3  | c^10
2  | c^15
```

You can easily access these elements using the following function:

```julia
springer_element(coxgroup(:E,7), 9)
```

which returns 

```julia
[1,4,6,2,3,5,7,1,4,6,2,3,5,7]
```

This is a reduced word for a $9$-Springer element of the finite Weyl group of type $E_7$.

The output of the function is intended to be fed directly into the previous three functions.

For example, if we are interested in the slope $\nu=4/9$ in type $E_7$, we would type:

```julia
count_points(coxgroup(:E,7), springer_element(coxgroup(:E,7), 9), 4)
```

```julia
interval_reps(coxgroup(:E,7), springer_element(coxgroup(:E,7), 9), 4)
```

```julia
count_points_lower(coxgroup(:E,7), springer_element(coxgroup(:E,7), 9), 4)
```

These replace hand-typing the long and error-prone function calls:

```julia
count_points(coxgroup(:E,7), [1,4,6,2,3,5,7,1,4,6,2,3,5,7], 4)
```

```julia
interval_reps(coxgroup(:E,7), [1,4,6,2,3,5,7,1,4,6,2,3,5,7], 4)
```

```julia
count_points_lower(coxgroup(:E,7), [1,4,6,2,3,5,7,1,4,6,2,3,5,7], 4)
```


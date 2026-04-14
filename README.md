# BraidStacks.jl

Tools to compute the number of points on braid stacks. This was written for the paper [arXiv:2603.20499](https://arxiv.org/abs/2603.20499). We heavily rely on Jean Michel's [Chevie](https://github.com/jmichel7/Chevie.jl). 


## Getting started

1. Download and install [Julia](https://julialang.org/downloads/). 
2. Install [Chevie](https://github.com/jmichel7/Chevie.jl):

```julia
using Pkg; Pkg.add("Chevie")
```

3. Copy-paste the script into Julia's command line.  
Alternatively, place `braid_stacks.jl` in Julia's bin folder then run:

```julia
include("braid_stacks.jl")
```


## Warm up example

Goal: Use `braid_stacks.jl` to determine the number of points on the braid stack $$M(β,γ)$$ for all $\gamma$ when $G = G_2$ and $β = (b_1b_2)^2$.  

We make the function call

```julia
count_points(coxgroup(:G,2),[1,2],2)
```

This prints the following:

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

We can study the isoclinic Deligne--Simpson problem by studying the point-counts of the braid stack $M(β,γ)$. Briefly, each slope $\nu$ gives rise to a braid $\beta_\nu$. Writing $\nu = d/m$ in lowest terms with $m$ a regular number for $W$, the braid looks like $$\beta_\nu = \widetilde{w}^d$$ for an $m$-Springer element $w\in W$, where $\widetilde{w}$ is the positive lift of $w$ to the braid group.

At the bottom of this readme, we give a detailed example of how to study the isoclinic Deligne--Simpson problem using the code. Beforehand, we give the relevant Weyl group elements $w$ for each regular number $m$ of each exceptional group $G$.  

When $m$ is regular elliptic, the element $w$ is of minimal length in its conjugacy class. We can use Chevie commands such as `classinfo(coxgroup(:G,2))` to find $w$, which has order $m$ and length $|\Phi|/m$.

When $m$ is regular but not elliptic, $w$ is determined in the appendix of:  

>Broué, Michel; Michel, Jean.  
Sur certains éléments réguliers des groupes de Weyl et les variétés de Deligne-Lusztig associées.(French)  
[Some regular elements of Weyl groups and the associated Deligne-Lusztig varieties]  
Finite reductive groups (Luminy, 1994), 73–139. Progr. Math., 141  

Here is the list:  

### $G = G_2$

| $m$ | $w$ |
| --- | --- |
| $6$ | $12 = c$ |
| $3$ | $c^2$ |
| $2$ | $c^3$ |

### $G = F_4$

| $m$ | $w$ |
| --- | --- |
| $12$ | $1234 = c$ |
| $8$ | $c23$ |
| $6$ | $c^2$ |
| $4$ | $c^3$ |
| $3$ | $c^4$ |
| $2$ | $c^6$ |

### $G = E_6$

| $m$ | $w$ |
| --- | --- |
| $12$ | $123456 = c$ |
| $9$ | $13432456$ |
| $8$ | $354163542 = x$ |
| $6$ | $c^2$ |
| $4$ | $x^2$ |
| $3$ | $c^4$ |
| $2$ | $x^4$ |

### $G = E_7$

| $m$ | $w$ |
| --- | --- |
| $18$ | $1234567 = c$ |
| $14$ | $42c^{-1}$ |
| $9$ | $(1462357)^2 = y^2$ |
| $7$ | $(134247654)^2$ |
| $6$ | $c^3$ |
| $3$ | $y^6$ |
| $2$ | $c^9$ |

### $G = E_8$

| $m$ | $w$ |
| --- | --- |
| $30$ | $12345678 = c$ |
| $24$ | $34c$ |
| $20$ | $4354c$ |
| $15$ | $c^2$ |
| $12$ | $2345c^2$ |
| $10$ | $c^3$ |
| $8$ | $123456c^3$ |
| $6$ | $c^5$ |
| $5$ | $c^6$ |
| $4$ | $(2314)^2 54234 56542 34576 54876 c^4$ |
| $3$ | $c^{10}$ |
| $2$ | $c^{15}$ |




## Examples using the table

Goal: Use `braid_stacks.jl` to compute the number of points on the braid stack $$M(\beta_\nu,\gamma)$$ for all unipotent classes $\gamma$ when $G = F_4$ and $\nu = 5/8$.  

In the list above, we have $w = c23 = 123423$ when $G = F_4$ and $m = 8$.  

This means $\beta_\nu = (c23)^5 = (123423)^5$. Therefore we make the function call

```julia
count_points(coxgroup(:F,4),[1,2,3,4,2,3],5)
```

(The last argument is $d=5$ because $\nu = d/m = 5/8$.)

This prints the following:

```julia
The group is G = F₄
The braid is β = b₁b₂b₁b₃b₂b₁b₃b₂b₃b₄b₃b₂b₁b₃b₂b₃b₄b₃b₂b₁b₃b₂b₃b₄b₁b₂b₃b₂b₄b₃
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

This means, in the notation of Jakob and Yun, we have $O_\nu = \widetilde{A_1}$ when $G=F_4$ and $\nu = 5/8$.  

To automatically find $O_\nu$, instead use 

```julia
interval_reps(coxgroup(:F,4),[1,2,3,4,2,3],5)
```

This prints the table above and returns the tuple

```julia
(UnipotentClass(Ã₁), UnipotentClass(F₄))
```

To count points at this lower representative $O_\nu = \widetilde{A_1}$, use

```julia
count_points_lower(coxgroup(:F,4),[1,2,3,4,2,3],5)
```

which returns the point-count in the row $ \gamma = O_\nu = \widetilde{A_1}$.




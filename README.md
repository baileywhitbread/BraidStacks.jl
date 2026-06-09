# BraidStacks.jl

Tools to compute the number of points on braid stacks. 

This was written for the paper [arXiv:2603.20499](https://arxiv.org/abs/2603.20499). 

We heavily rely on Jean Michel's [Chevie](https://github.com/jmichel7/Chevie.jl). 


## Getting started

1. Download and install [Julia](https://julialang.org/downloads/). 
2. Install [Chevie](https://github.com/jmichel7/Chevie.jl):

```julia
using Pkg; Pkg.add("Chevie")
```

Make sure you are using Chevie v0.1.14 or greater. 

Check:

```julia
using Pkg; Pkg.status("Chevie")
```

Upgrade:

```julia
using Pkg; Pkg.add(url="https://github.com/jmichel7/Chevie.jl.git#v0.1.14")
```



3. Copy-paste the script into Julia's command line.  
Alternatively, place `braid_stacks.jl` in Julia's bin folder then run:

```julia
include("braid_stacks.jl")
```


## Warm up example

Goal: Count points on $$M(β,γ)$$ for all unipotent classes $\gamma$ when $G = G_2$ and $β = (b_1b_2)^2$.  

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

We solve the isoclinic Deligne--Simpson problem by counting points on the braid stack $M(β,γ)$. Each slope $\nu$ gives rise to a braid $\beta_\nu$. Write $\nu = d/m$ in lowest terms with $m$ a regular number for $W$. Up to cyclic shift, the braid looks like $$\beta_\nu = \widetilde{w}^d$$ for an $m$-Springer element $w\in W$, where $\widetilde{w}$ is the positive lift of $w$ to the braid group.

Goal: Count points on $$M(\beta_\nu,\gamma)$$ for all unipotent classes $\gamma$ when $G = F_4$ and $\nu = 5/8$.  

In next section, you will find that $w = 142323$ is an $m$-Springer element when $G = F_4$ and $m = 8$.  

This means $\beta_\nu = (142323)^5$. So we call the function

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

To automatically find the minimal unipotent class (if it exists!), instead use 

```julia
interval_reps(coxgroup(:F,4),[1,4,2,3,2,3],5)
```

This prints the table above and returns the tuple

```julia
(UnipotentClass(Ã₁), UnipotentClass(F₄))
```

To count points at this lower representative $C_\nu = \widetilde{A_1}$, use

```julia
count_points_lower(coxgroup(:F,4),[1,4,2,3,2,3],5)
```

which returns the point-count in the row $\gamma = C_\nu = \widetilde{A_1}$.





## Tables of Springer elements

For reference, we give a table of $m$-Springer elements $w$ for each regular number $m$ of each exceptional group $G$, using the appendix of:  

>Broué, Michel; Michel, Jean.  
Sur certains éléments réguliers des groupes de Weyl et les variétés de Deligne-Lusztig associées.(French)  
[Some regular elements of Weyl groups and the associated Deligne-Lusztig varieties]  
Finite reductive groups (Luminy, 1994), 73–139. Progr. Math., 141  

Here is the list:  

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
[1,3,4,2,4,7,6,5,4,1,3,4,2,4,7,6,5,4]
```


The output of this function is intended to be fed directly into the previous three functions, for example as follows:

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
count_points(coxgroup(:E,7), [1,3,4,2,4,7,6,5,4,1,3,4,2,4,7,6,5,4], 4)
```

```julia
interval_reps(coxgroup(:E,7), [1,3,4,2,4,7,6,5,4,1,3,4,2,4,7,6,5,4], 4)
```

```julia
count_points_lower(coxgroup(:E,7), [1,3,4,2,4,7,6,5,4,1,3,4,2,4,7,6,5,4], 4)
```

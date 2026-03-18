# BraidStacks.jl

Tools to compute the number of points on braid stacks. We heavily rely on Jean Michel's [Chevie](https://github.com/jmichel7/Chevie.jl). 


## Getting started

1. Download and install [Julia](https://julialang.org/downloads/). 
2. Install [Chevie](https://github.com/jmichel7/Chevie.jl):

```julia
using Pkg; Pkg.add("Chevie")
```
3. Place the braid_stacks.jl file in Julia's bin folder. In Julia's command-line, run the below:

```julia
include("braid_stacks.jl")
```

3. (Alternative to step 3 above). Copy-paste the script into Julia's command line.


## Examples

Goal: Use `braid_stacks.jl` to determine the number of points on the braid stack $M(β,γ)$ when $G = G_2$ and $β = (b_1b_2)^2$.  

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






## A list of Weyl group elements

To solve the isoclinic Deligne--Simpson problem, one studies the braid stack $M(β,γ)$. To do so, each slope $\nu$ gives rise to a braid $\beta_\nu$. As mentioned in the preprint, writing $\nu d/m$ in lowest terms with $m$ a regular number for $W$, the braid looks like $\widetilde{w}^d$ for an $m$-Springer element $w\in W$. (Here, $\widetilde{w}$ denotes the lift of $w$ to the positive braid $\widetilde{w}$ in the positive braid monoid.)

Below we give lists of such w for each regular number $m$ of each exceptional group $G$. When $m$ is elliptic, it will be of minimal length in its conjugacy class. Therefore we can use Chevie commands such as `classinfo(coxgroup(:G,2))` to find $w$. Specifically, $w$ has order $m$ and length $|\Phi|/m$.

When $m$ is not elliptic, $w$ is determined in the appendix of:  

>Broué, Michel; Michel, Jean.  
Sur certains éléments réguliers des groupes de Weyl et les variétés de Deligne-Lusztig associées.(French)  
[Some regular elements of Weyl groups and the associated Deligne-Lusztig varieties]  
Finite reductive groups (Luminy, 1994), 73–139. Progr. Math., 141  

Here is the list:

### G = G2

| m | w |
| --- | --- |
| 6 | `12 = c` |
| 3 | `c^2` |
| 2 | `c^3` |

### G = F4

| m | w |
| --- | --- |
| 12 | `1234 = c` |
| 8 | `c23` |
| 6 | `c^2` |
| 4 | `c^3` |
| 3 | `c^4` |
| 2 | `c^6` |

### G = E6

| m | w |
| --- | --- |
| 12 | `123456 = c` |
| 9 | `13432456` |
| 8 | `354163542 = x` |
| 6 | `c^2` |
| 4 | `x^2` |
| 3 | `c^4` |
| 2 | `x^4` |

### G = E7

| m | w |
| --- | --- |
| 18 | `1234567 = c` |
| 14 | `42c^(-1)` |
| 9 | `(1462357)^2 = y^2` |
| 7 | `(134247654)^2` |
| 6 | `c^3` |
| 3 | `y^6` |
| 2 | `c^9` |

### G = E8

| m | w |
| --- | --- |
| 30 | `12345678 = c` |
| 24 | `34c` |
| 20 | `4354c` |
| 15 | `c^2` |
| 12 | `2345c^2` |
| 10 | `c^3` |
| 8 | `123456c^3` |
| 6 | `c^5` |
| 5 | `c^6` |
| 4 | `(2314)^2 54234 56542 34576 54876 c^4` |
| 3 | `c^10` |
| 2 | `c^15` |



## Examples using the table

Goal: Use `braid_stacks.jl` to determine $O_\nu$ for $G = F_4$ and $\nu = 5/8$.  

In the list above, we have $w = c23 = 123423$ when $G = F_4$ and $m = 8$. 

Therefore we make the function call

```julia
count_points(coxgroup(:F,4),[1,2,3,4,2,3],5)
```

One visually inspects the table and finds $O_\nu = \widetilde{A_1}$.  

(This is often difficult to determine, unless you have the poset memorised!)  

To automatically find $O_\nu$, instead use 

```julia
interval_reps(coxgroup(:F,4),[1,2,3,4,2,3],5)
```


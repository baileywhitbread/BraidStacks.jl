# BraidStacks.jl

Tools to compute the number of points on braid stacks and address the isoclinic Deligne--Simpson problem. 
This was written for the paper [arXiv:2603.20499](https://arxiv.org/abs/2603.20499). 
We heavily rely on Jean Michel's port of the computer algebra system [Chevie.jl](https://github.com/jmichel7/Chevie.jl).

![](https://github.com/baileywhitbread/BraidStacks.jl/blob/main/BraidStacksAnimation.gif)

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

3. Paste the contents of `braid_stacks_stable.jl` into Julia's command line.  
Alternatively, place the file `braid_stacks_stable.jl` in Julia's working directory, then run:

```julia
include("braid_stacks_stable.jl")
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

## Proofs regarding the isoclinic Deligne--Simpson problem

The script `braid_stacks_proofs.jl` contains proofs for theorems in:
>Kamgarpour, Masoud; Whitbread, Bailey.  
The isoclinic Deligne--Simpson problem and rigid connections
https://arxiv.org/pdf/2603.20499 

Namely, the script proves:
1. The correctness of Table 1 in Theorem 1.
2. The physical rigidity of the connections in Theorem 2 (a).

More specifically, running the script prints 

```julia
G = G₂, ν = 1//2: C_ν = UnipotentClass(Ã₁)
G = G₂, ν = 1//3: C_ν = UnipotentClass(G₂(a₁))
G = G₂, ν = 2//3: C_ν = UnipotentClass(A₁)
G = G₂, ν = 1//6: C_ν = UnipotentClass(G₂)
G = G₂, ν = 5//6: C_ν = UnipotentClass(A₁)
G = F₄, ν = 1//2: C_ν = UnipotentClass(A₁+Ã₁)
G = F₄, ν = 1//3: C_ν = UnipotentClass(Ã₂+A₁)
G = F₄, ν = 2//3: C_ν = UnipotentClass(Ã₁)
G = F₄, ν = 1//4: C_ν = UnipotentClass(F₄(a₃))
G = F₄, ν = 3//4: C_ν = UnipotentClass(A₁)
G = F₄, ν = 1//6: C_ν = UnipotentClass(F₄(a₂))
G = F₄, ν = 5//6: C_ν = UnipotentClass(A₁)
G = F₄, ν = 1//8: C_ν = UnipotentClass(F₄(a₁))
G = F₄, ν = 3//8: C_ν = UnipotentClass(A₂+Ã₁)
G = F₄, ν = 5//8: C_ν = UnipotentClass(Ã₁)
G = F₄, ν = 7//8: C_ν = UnipotentClass(A₁)
G = F₄, ν = 1//12: C_ν = UnipotentClass(F₄)
G = F₄, ν = 5//12: C_ν = UnipotentClass(A₂+Ã₁)
G = F₄, ν = 7//12: C_ν = UnipotentClass(A₁+Ã₁)
G = F₄, ν = 11//12: C_ν = UnipotentClass(A₁)
G = E₆, ν = 1//2: C_ν = UnipotentClass(3A₁)
G = E₆, ν = 1//3: C_ν = UnipotentClass(2A₂+A₁)
G = E₆, ν = 2//3: C_ν = UnipotentClass(2A₁)
G = E₆, ν = 1//4: C_ν = UnipotentClass(D₄(a₁))
G = E₆, ν = 3//4: C_ν = UnipotentClass(A₁)
G = E₆, ν = 1//6: C_ν = UnipotentClass(E₆(a₃))
G = E₆, ν = 5//6: C_ν = UnipotentClass(A₁)
G = E₆, ν = 1//8: C_ν = UnipotentClass(D₅)
G = E₆, ν = 3//8: C_ν = UnipotentClass(A₂+2A₁)
G = E₆, ν = 5//8: C_ν = UnipotentClass(2A₁)
G = E₆, ν = 7//8: C_ν = UnipotentClass(A₁)
G = E₆, ν = 1//9: C_ν = UnipotentClass(E₆(a₁))
G = E₆, ν = 2//9: C_ν = UnipotentClass(A₄+A₁)
G = E₆, ν = 4//9: C_ν = UnipotentClass(A₂+A₁)
G = E₆, ν = 5//9: C_ν = UnipotentClass(3A₁)
G = E₆, ν = 7//9: C_ν = UnipotentClass(A₁)
G = E₆, ν = 8//9: C_ν = UnipotentClass(A₁)
G = E₆, ν = 1//12: C_ν = UnipotentClass(E₆)
G = E₆, ν = 5//12: C_ν = UnipotentClass(A₂+2A₁)
G = E₆, ν = 7//12: C_ν = UnipotentClass(3A₁)
G = E₆, ν = 11//12: C_ν = UnipotentClass(A₁)
G = E₇, ν = 1//2: C_ν = UnipotentClass(4A₁)
G = E₇, ν = 1//3: C_ν = UnipotentClass(2A₂+A₁)
G = E₇, ν = 2//3: C_ν = UnipotentClass(2A₁)
G = E₇, ν = 1//6: C_ν = UnipotentClass(E₇(a₅))
G = E₇, ν = 5//6: C_ν = UnipotentClass(A₁)
G = E₇, ν = 1//7: C_ν = UnipotentClass(A₆)
G = E₇, ν = 2//7: C_ν = UnipotentClass(A₃+A₂)
G = E₇, ν = 3//7: C_ν = UnipotentClass(A₂+2A₁)
G = E₇, ν = 4//7: C_ν = UnipotentClass(3A₁′)
G = E₇, ν = 5//7: C_ν = UnipotentClass(2A₁)
G = E₇, ν = 6//7: C_ν = UnipotentClass(A₁)
G = E₇, ν = 1//9: C_ν = UnipotentClass(E₆(a₁))
G = E₇, ν = 2//9: C_ν = UnipotentClass(A₄+A₁)
G = E₇, ν = 4//9: C_ν = UnipotentClass(A₂+A₁)
G = E₇, ν = 5//9: C_ν = UnipotentClass(3A₁′)
G = E₇, ν = 7//9: C_ν = UnipotentClass(A₁)
G = E₇, ν = 8//9: C_ν = UnipotentClass(A₁)
G = E₇, ν = 1//14: C_ν = UnipotentClass(E₇(a₁))
G = E₇, ν = 3//14: C_ν = UnipotentClass(A₄+A₂)
G = E₇, ν = 5//14: C_ν = UnipotentClass(2A₂+A₁)
G = E₇, ν = 9//14: C_ν = UnipotentClass(2A₁)
G = E₇, ν = 11//14: C_ν = UnipotentClass(A₁)
G = E₇, ν = 13//14: C_ν = UnipotentClass(A₁)
G = E₇, ν = 1//18: C_ν = UnipotentClass(E₇)
G = E₇, ν = 5//18: C_ν = UnipotentClass(A₃+A₂+A₁)
G = E₇, ν = 7//18: C_ν = UnipotentClass(A₂+3A₁)
G = E₇, ν = 11//18: C_ν = UnipotentClass(3A₁′)
G = E₇, ν = 13//18: C_ν = UnipotentClass(2A₁)
G = E₇, ν = 17//18: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//2: C_ν = UnipotentClass(4A₁)
G = E₈, ν = 1//3: C_ν = UnipotentClass(2A₂+2A₁)
G = E₈, ν = 2//3: C_ν = UnipotentClass(2A₁)
G = E₈, ν = 1//4: C_ν = UnipotentClass(2A₃)
G = E₈, ν = 3//4: C_ν = UnipotentClass(2A₁)
G = E₈, ν = 1//5: C_ν = UnipotentClass(A₄+A₃)
G = E₈, ν = 2//5: C_ν = UnipotentClass(A₂+3A₁)
G = E₈, ν = 3//5: C_ν = UnipotentClass(3A₁)
G = E₈, ν = 4//5: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//6: C_ν = UnipotentClass(E₈(a₇))
G = E₈, ν = 5//6: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//8: C_ν = UnipotentClass(A₇)
G = E₈, ν = 3//8: C_ν = UnipotentClass(2A₂+A₁)
G = E₈, ν = 5//8: C_ν = UnipotentClass(3A₁)
G = E₈, ν = 7//8: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//10: C_ν = UnipotentClass(E₈(a₆))
G = E₈, ν = 3//10: C_ν = UnipotentClass(D₄(a₁)+A₁)
G = E₈, ν = 7//10: C_ν = UnipotentClass(2A₁)
G = E₈, ν = 9//10: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//12: C_ν = UnipotentClass(E₈(a₅))
G = E₈, ν = 5//12: C_ν = UnipotentClass(A₂+3A₁)
G = E₈, ν = 7//12: C_ν = UnipotentClass(3A₁)
G = E₈, ν = 11//12: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//15: C_ν = UnipotentClass(E₈(a₄))
G = E₈, ν = 2//15: C_ν = UnipotentClass(D₇(a₂))
G = E₈, ν = 4//15: C_ν = UnipotentClass(D₄(a₁)+A₂)
G = E₈, ν = 7//15: C_ν = UnipotentClass(A₂+A₁)
G = E₈, ν = 8//15: C_ν = UnipotentClass(4A₁)
G = E₈, ν = 11//15: C_ν = UnipotentClass(2A₁)
G = E₈, ν = 13//15: C_ν = UnipotentClass(A₁)
G = E₈, ν = 14//15: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//20: C_ν = UnipotentClass(E₈(a₂))
G = E₈, ν = 3//20: C_ν = UnipotentClass(A₆+A₁)
G = E₈, ν = 7//20: C_ν = UnipotentClass(2A₂+2A₁)
G = E₈, ν = 9//20: C_ν = UnipotentClass(A₂+2A₁)
G = E₈, ν = 11//20: C_ν = UnipotentClass(4A₁)
G = E₈, ν = 13//20: C_ν = UnipotentClass(2A₁)
G = E₈, ν = 17//20: C_ν = UnipotentClass(A₁)
G = E₈, ν = 19//20: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//24: C_ν = UnipotentClass(E₈(a₁))
G = E₈, ν = 5//24: C_ν = UnipotentClass(A₄+A₃)
G = E₈, ν = 7//24: C_ν = UnipotentClass(A₃+A₂+A₁)
G = E₈, ν = 11//24: C_ν = UnipotentClass(A₂+2A₁)
G = E₈, ν = 13//24: C_ν = UnipotentClass(4A₁)
G = E₈, ν = 17//24: C_ν = UnipotentClass(2A₁)
G = E₈, ν = 19//24: C_ν = UnipotentClass(A₁)
G = E₈, ν = 23//24: C_ν = UnipotentClass(A₁)
G = E₈, ν = 1//30: C_ν = UnipotentClass(E₈)
G = E₈, ν = 7//30: C_ν = UnipotentClass(A₄+A₂+A₁)
G = E₈, ν = 11//30: C_ν = UnipotentClass(2A₂+2A₁)
G = E₈, ν = 13//30: C_ν = UnipotentClass(A₂+3A₁)
G = E₈, ν = 17//30: C_ν = UnipotentClass(4A₁)
G = E₈, ν = 19//30: C_ν = UnipotentClass(3A₁)
G = E₈, ν = 23//30: C_ν = UnipotentClass(2A₁)
G = E₈, ν = 29//30: C_ν = UnipotentClass(A₁)
```

and

```julia
G = A₁, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = A₁, ν = 3//2: |M(β_ν,C_ν)^F| = 1
G = A₂, ν = 1//3: |M(β_ν,C_ν)^F| = 1
G = A₂, ν = 4//3: |M(β_ν,C_ν)^F| = 1
G = A₂, ν = 2//3: |M(β_ν,C_ν)^F| = 1
G = A₃, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = A₃, ν = 5//4: |M(β_ν,C_ν)^F| = 1
G = A₃, ν = 3//4: |M(β_ν,C_ν)^F| = 1
G = A₄, ν = 1//5: |M(β_ν,C_ν)^F| = 1
G = A₄, ν = 6//5: |M(β_ν,C_ν)^F| = 1
G = A₄, ν = 2//5: |M(β_ν,C_ν)^F| = 1
G = A₄, ν = 3//5: |M(β_ν,C_ν)^F| = 1
G = A₄, ν = 4//5: |M(β_ν,C_ν)^F| = 1
G = A₅, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = A₅, ν = 7//6: |M(β_ν,C_ν)^F| = 1
G = A₅, ν = 5//6: |M(β_ν,C_ν)^F| = 1
G = A₆, ν = 1//7: |M(β_ν,C_ν)^F| = 1
G = A₆, ν = 8//7: |M(β_ν,C_ν)^F| = 1
G = A₆, ν = 2//7: |M(β_ν,C_ν)^F| = 1
G = A₆, ν = 3//7: |M(β_ν,C_ν)^F| = 1
G = A₆, ν = 4//7: |M(β_ν,C_ν)^F| = 1
G = A₆, ν = 6//7: |M(β_ν,C_ν)^F| = 1
G = A₇, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = A₇, ν = 9//8: |M(β_ν,C_ν)^F| = 1
G = A₇, ν = 3//8: |M(β_ν,C_ν)^F| = 1
G = A₇, ν = 7//8: |M(β_ν,C_ν)^F| = 1
G = A₈, ν = 1//9: |M(β_ν,C_ν)^F| = 1
G = A₈, ν = 10//9: |M(β_ν,C_ν)^F| = 1
G = A₈, ν = 2//9: |M(β_ν,C_ν)^F| = 1
G = A₈, ν = 4//9: |M(β_ν,C_ν)^F| = 1
G = A₈, ν = 5//9: |M(β_ν,C_ν)^F| = 1
G = A₈, ν = 8//9: |M(β_ν,C_ν)^F| = 1
G = B₂, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = B₂, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = B₂, ν = 5//4: |M(β_ν,C_ν)^F| = 1
G = B₂, ν = 3//4: |M(β_ν,C_ν)^F| = 1
G = B₃, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = B₃, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = B₃, ν = 7//6: |M(β_ν,C_ν)^F| = 1
G = B₄, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = B₄, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = B₄, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = B₄, ν = 9//8: |M(β_ν,C_ν)^F| = 1
G = B₄, ν = 3//8: |M(β_ν,C_ν)^F| = 1
G = B₄, ν = 5//8: |M(β_ν,C_ν)^F| = 1
G = B₄, ν = 3//4: |M(β_ν,C_ν)^F| = 1
G = B₅, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = B₅, ν = 1//10: |M(β_ν,C_ν)^F| = 1
G = B₅, ν = 11//10: |M(β_ν,C_ν)^F| = 1
G = B₅, ν = 3//10: |M(β_ν,C_ν)^F| = 1
G = B₆, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = B₆, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = B₆, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = B₆, ν = 1//12: |M(β_ν,C_ν)^F| = 1
G = B₆, ν = 13//12: |M(β_ν,C_ν)^F| = 1
G = B₆, ν = 7//12: |M(β_ν,C_ν)^F| = 1
G = B₇, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = B₇, ν = 1//14: |M(β_ν,C_ν)^F| = 1
G = B₇, ν = 15//14: |M(β_ν,C_ν)^F| = 1
G = B₇, ν = 3//14: |M(β_ν,C_ν)^F| = 1
G = B₇, ν = 5//14: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 1//16: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 17//16: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 3//16: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 9//16: |M(β_ν,C_ν)^F| = 1
G = B₈, ν = 3//8: |M(β_ν,C_ν)^F| = 1
G = C₃, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = C₃, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = C₃, ν = 7//6: |M(β_ν,C_ν)^F| = 1
G = C₃, ν = 5//6: |M(β_ν,C_ν)^F| = 1
G = C₄, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = C₄, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = C₄, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = C₄, ν = 9//8: |M(β_ν,C_ν)^F| = 1
G = C₄, ν = 3//8: |M(β_ν,C_ν)^F| = 1
G = C₄, ν = 7//8: |M(β_ν,C_ν)^F| = 1
G = C₅, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = C₅, ν = 1//10: |M(β_ν,C_ν)^F| = 1
G = C₅, ν = 11//10: |M(β_ν,C_ν)^F| = 1
G = C₅, ν = 3//10: |M(β_ν,C_ν)^F| = 1
G = C₅, ν = 9//10: |M(β_ν,C_ν)^F| = 1
G = C₆, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = C₆, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = C₆, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = C₆, ν = 1//12: |M(β_ν,C_ν)^F| = 1
G = C₆, ν = 13//12: |M(β_ν,C_ν)^F| = 1
G = C₆, ν = 11//12: |M(β_ν,C_ν)^F| = 1
G = C₇, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = C₇, ν = 1//14: |M(β_ν,C_ν)^F| = 1
G = C₇, ν = 15//14: |M(β_ν,C_ν)^F| = 1
G = C₇, ν = 3//14: |M(β_ν,C_ν)^F| = 1
G = C₇, ν = 5//14: |M(β_ν,C_ν)^F| = 1
G = C₇, ν = 13//14: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 1//16: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 17//16: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 3//16: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 5//16: |M(β_ν,C_ν)^F| = 1
G = C₈, ν = 15//16: |M(β_ν,C_ν)^F| = 1
G = D₄, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = D₄, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = D₄, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = D₄, ν = 7//6: |M(β_ν,C_ν)^F| = 1
G = D₄, ν = 3//4: |M(β_ν,C_ν)^F| = 1
G = D₅, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = D₅, ν = 9//8: |M(β_ν,C_ν)^F| = 1
G = D₅, ν = 3//8: |M(β_ν,C_ν)^F| = 1
G = D₅, ν = 5//8: |M(β_ν,C_ν)^F| = 1
G = D₆, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = D₆, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = D₆, ν = 1//10: |M(β_ν,C_ν)^F| = 1
G = D₆, ν = 11//10: |M(β_ν,C_ν)^F| = 1
G = D₆, ν = 3//10: |M(β_ν,C_ν)^F| = 1
G = D₇, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = D₇, ν = 1//12: |M(β_ν,C_ν)^F| = 1
G = D₇, ν = 13//12: |M(β_ν,C_ν)^F| = 1
G = D₇, ν = 7//12: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 1//14: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 15//14: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 3//14: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 5//14: |M(β_ν,C_ν)^F| = 1
G = D₈, ν = 3//8: |M(β_ν,C_ν)^F| = 1
G = G₂, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = G₂, ν = 1//3: |M(β_ν,C_ν)^F| = 1
G = G₂, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = G₂, ν = 7//6: |M(β_ν,C_ν)^F| = 1
G = G₂, ν = 2//3: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 1//3: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 1//12: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 13//12: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 3//8: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 5//8: |M(β_ν,C_ν)^F| = 1
G = F₄, ν = 3//4: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 1//3: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 1//9: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 1//12: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 13//12: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 2//9: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 4//9: |M(β_ν,C_ν)^F| = 1
G = E₆, ν = 7//9: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 1//14: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 1//18: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 19//18: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 3//14: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 9//14: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 11//14: |M(β_ν,C_ν)^F| = 1
G = E₇, ν = 7//18: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//2: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//3: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//4: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//5: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//6: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//8: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//10: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//12: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//15: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//20: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//24: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 1//30: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 31//30: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 3//10: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 2//15: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 4//15: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 7//15: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 3//20: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 13//20: |M(β_ν,C_ν)^F| = 1
G = E₈, ν = 19//24: |M(β_ν,C_ν)^F| = 1
```
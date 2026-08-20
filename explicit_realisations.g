# Load Willem de Graaf's package for simple Lie algebras
# https://gap-packages.github.io/sla/
# and load the exceptional simple Lie algebras.

LoadPackage("sla");;
G2 := SimpleLieAlgebra("G",2,Rationals);;
F4 := SimpleLieAlgebra("F",4,Rationals);;
E6 := SimpleLieAlgebra("E",6,Rationals);;
E7 := SimpleLieAlgebra("E",7,Rationals);;
E8 := SimpleLieAlgebra("E",8,Rationals);;

# Below are the pairs [x,D] where
# x = 1,2,3,...,h-1 is coprime to the Coxeter number h
# D is the weighted Dynkin diagram for C_{x/h} given in Appendix A of

# de Graaf, Willem A.
# Computing with nilpotent orbits in simple Lie algebras of exceptional type.
# LMS J. Comput. Math. 11 (2008), 280–297.

# Beside each pair [x,D] is a comment matching D with C_{x/h} according to de Graaf's Appendix A.
# (One can use DisplayWeightedDynkinDiagram for this.)

cases := [
  [G2,[
  [1,[2,2]], 				# [2,2] = G2 = C_{1/6}
  [5,[0,1]]					# [0,1] = A1 = C_{5/6}
  ]],
  
  [F4,[
  [1,[2,2,2,2]], 			# [2,2,2,2] = F4 			= C_{1/12}
  [5,[0,0,1,0]], 			# [0,0,1,0] = A1' + A2 		= C_{5/12}
  [7,[0,0,0,1]], 			# [0,0,0,1] = A1 + A1' 		= C_{7/12}
  [11,[0,1,0,0]]			# [0,1,0,0] = A1 			= C_{11/12}
  ]],
  
  [E6,[
  [1,[2,2,2,2,2,2]], 		# [2,2,2,2,2,2] = E6 			= C_{1/12}
  [5,[0,0,1,0,1,0]], 		# [0,0,1,0,1,0] = 2A1 + A2 		= C_{5/12}
  [7,[0,0,0,1,0,0]], 		# [0,0,0,1,0,0] = 3A1 			= C_{7/12}
  [11,[0,1,0,0,0,0]]		# [0,1,0,0,0,0] = A1 			= C_{11/12}
  ]],
  
  [E7,[
  [1,[2,2,2,2,2,2,2]], 		# [2,2,2,2,2,2,2] = E7 				= C_{1/18}
  [5,[0,0,0,0,2,0,0]], 		# [0,0,0,0,2,0,0] = A3 + A2 + A1 	= C_{5/18}
  [7,[0,2,0,0,0,0,0]], 		# [0,2,0,0,0,0,0] = A2 + 3A1 		= C_{7/18}
  [11,[0,0,1,0,0,0,0]], 	# [0,0,1,0,0,0,0] = (3A1)' 			= C_{11/18}
  [13,[0,0,0,0,0,1,0]], 	# [0,0,0,0,0,1,0] = 2A1 			= C_{13/18}
  [17,[1,0,0,0,0,0,0]]		# [1,0,0,0,0,0,0] = A1 				= C_{17/18}
  ]],
  
  [E8,[
  [1,[2,2,2,2,2,2,2,2]], 	# [2,2,2,2,2,2,2,2] = E8				= C_{1/30}
  [7,[0,0,1,0,0,1,0,0]], 	# [0,0,1,0,0,1,0,0] = A4 + A2 + A1		= C_{7/30}
  [11,[0,0,0,0,1,0,0,0]], 	# [0,0,0,0,1,0,0,0] = 2A2 + 2A1			= C_{11/30}
  [13,[0,0,1,0,0,0,0,0]], 	# [0,0,1,0,0,0,0,0] = A2 + 3A1			= C_{13/30}
  [17,[0,1,0,0,0,0,0,0]], 	# [0,1,0,0,0,0,0,0] = 4A1				= C_{17/30}
  [19,[0,0,0,0,0,0,1,0]], 	# [0,0,0,0,0,0,1,0] = 3A1				= C_{19/30}
  [23,[1,0,0,0,0,0,0,0]], 	# [1,0,0,0,0,0,0,0] = 2A1				= C_{23/30}
  [29,[0,0,0,0,0,0,0,1]]	# [0,0,0,0,0,0,0,1] = A1				= C_{29/30}
  ]]
];;

# Now let us check that the weighted Dynkin diagram of N_x is D
# for each exceptional simple Lie algebra and
# for each x = 1,2,3,...,h-1 coprime to h.

for case in cases do
  L := case[1];;
  R := RootSystem(L);;
  PosR := PositiveRootsNF(R);;
  v := NegativeRootVectors(R);;

  for x_D_pair in case[2] do
	x := x_D_pair[1];;
	D := x_D_pair[2];;
    I := Filtered([1..Length(PosR)],i->Sum(PosR[i])=x);;
    N_x := Sum(v{I});;
    Assert(0,WeightedDynkinDiagram(L,N_x)=D);
  od;
od;

Print("All cases passed.\n");

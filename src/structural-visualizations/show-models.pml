
load /Users/idamei/polyphenol-oxidases/data/alphafold/AF-A0A1C9CXH9-F1-model_v4-MtPPO7.pdb, XP_003662515.1
load /Users/idamei/polyphenol-oxidases/data/alphafold/AF-G2Q526-F1-model_v4-MtPPO-809.pdb, XP_003659809.1
load /Users/idamei/polyphenol-oxidases/data/alphafold/AF-Q2GZJ4-F1-model_v4_CgPPO-266.pdb, XP_001224266.1
load /Users/idamei/polyphenol-oxidases/data/alphafold/AF-Q2H7I7-F1-model_v4_CgPPO-1473.pdb, XP_001221473.1
load /Users/idamei/polyphenol-oxidases/data/alphafold/c2092.pdb, c2092_g1_i1
fetch 6Z1S, XP_003666010.1
fetch 2y9x
fetch 2y9w
fetch 4j3p

select chain_D_2y9x, (chain D and 2y9x)
hide cartoon, (!chain_D_2y9x and 2y9x)
hide spheres, (!chain_D_2y9x and 2y9x)
hide sticks, (!chain_D_2y9x and 2y9x)

select chain_B_2y9w, (chain B and 2y9w)
hide cartoon, (!chain_B_2y9w and 2y9w)
hide spheres, (!chain_B_2y9w and 2y9w)
hide sticks, (!chain_B_2y9w and 2y9w) 

cealign XP_003662515.1, XP_003659809.1
cealign XP_003662515.1, XP_001224266.1
cealign XP_003662515.1, XP_001221473.1
cealign XP_003662515.1, c2092_g1_i1
cealign XP_003662515.1, XP_003666010.1
cealign XP_003662515.1, 2y9x
cealign XP_003662515.1, 2y9w
cealign XP_003662515.1, 4j3p

remove solvent

select cobber_1, (chain_B_2y9w and resi 401)
select cobber_2, (chain_B_2y9w and resi 400)

select as_2y9w, 2y9w within 9 of cobber_1
show sticks, as_2y9w
label n.  CA and as_2y9w, "%s-%s" % (resn,resi)
select as_2y9x, 2y9x within 9 of cobber_1
show sticks, as_2y9x
label n.  CA and as_2y9x, "%s-%s" % (resn,resi)
select as_XP_003662515.1, (XP_003662515.1 within 9 of cobber_1)
show sticks, as_XP_003662515.1
label n.  CA and as_XP_003662515.1, "%s-%s" % (resn,resi)
select as_XP_003659809.1, (XP_003659809.1 within 9 of cobber_1)
show sticks, as_XP_003659809.1
label n.  CA and as_XP_003659809.1, "%s-%s" % (resn,resi)
select as_XP_001224266.1, XP_001224266.1 within 9 of cobber_1
show sticks, as_XP_001224266.1
label n.  CA and as_XP_001224266.1, "%s-%s" % (resn,resi)
select as_XP_001221473.1, XP_001221473.1 within 9 of cobber_1
show sticks, as_XP_001221473.1
label n.  CA and as_XP_001221473.1, "%s-%s" % (resn,resi)
select as_c2092_g1_i1, c2092_g1_i1 within 9 of cobber_1
show sticks, as_c2092_g1_i1
label n.  CA and as_c2092_g1_i1, "%s-%s" % (resn,resi)
select as_XP_003666010.1, XP_003666010.1 within 9 of cobber_1
show sticks, as_XP_003666010.1
label n.  CA and as_XP_003666010.1, "%s-%s" % (resn,resi)
select as_4j3p, 4j3p within 9 of cobber_1
show sticks, as_4j3p
label n.  CA and as_4j3p, "%s-%s" % (resn,resi)

@src/structural-visualizations/nicify.pml
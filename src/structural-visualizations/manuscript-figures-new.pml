@src/pymol-visualizations/nicify.pml

load data/compare-architectures/pdbs/4j3p.pdb, 4j3p
color orange, 4j3p
select cons_4j3p, resi 102 and 4j3p or resi 110 and 4j3p or resi 119 and 4j3p or resi 239 and 4j3p or resi 284 and 4j3p or resi 288 and 4j3p or resi 311 and 4j3p or resi 273 and 4j3p or resi 285 and 4j3p
show licorice, cons_4j3p
color atomic, (cons_4j3p and not elem C)
label n. CA and cons_4j3p, "%s-%s" % (resn, resi) 

load data/compare-architectures/pdbs/2y9x.pdb, 2y9x
color teal, 2y9x
select cons_2y9x, resi 61 and 2y9x or resi 85 and 2y9x or resi 94 and 2y9x or resi 200 and 2y9x or resi 256 and 2y9x or resi 259 and 2y9x or resi 263 and 2y9x or resi 295 and 2y9x or resi 260 and 2y9x
show licorice, cons_2y9x
color atomic, (cons_2y9x not elem C)
label n. CA and cons_2y9x, "%s-%s" % (resn, resi) 
select tropolone, 2y9x and resi 410
color atomic, tropolone
color yellow, (tropolone and elem C)

remove chain B or chain C or chain D or chain E or chain F or chain G or chain H or chain I or chain J

cealign 2y9x, 4j3p

hide (resi 2011 and 4j3p)

set sphere_scale, 0.5

@src/structural-visualizations/nicify.pml
set cartoon_transparency, 0.8
hide labels

set_view (\
     0.343452483,   -0.027162792,   -0.938776255,\
     0.225473017,   -0.967959523,    0.110496975,\
    -0.911700845,   -0.249622017,   -0.326327145,\
    -0.000027796,   -0.000036873,  -73.859146118,\
    -5.644330978,  -26.331651688,  -37.025714874,\
    33.452472687,  114.269111633,  -20.000000000 )

disable
enable 4j3p
ray
png manuscript/figures/4j3p.png

disable
enable 2y9x
ray
png manuscript/figures/2y9x.png
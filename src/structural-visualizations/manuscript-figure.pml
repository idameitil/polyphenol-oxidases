@src/pymol-visualizations/nicify.pml

load data/compare-architectures/pdbs/4j3p.pdb, 4j3p
color orange, 4j3p
select cons_4j3p, resi 102 and 4j3p or resi 110 and 4j3p or resi 119 and 4j3p or resi 239 and 4j3p or resi 284 and 4j3p or resi 288 and 4j3p or resi 311 and 4j3p or resi 302 and 4j3p
show licorice, cons_4j3p
color atomic, (cons_4j3p and not elem C)
label n. CA and cons_4j3p, "%s-%s" % (resn, resi) 

load data/compare-architectures/pdbs/2y9x.pdb, 2y9x
color teal, 2y9x
select cons_2y9x, resi 61 and 2y9x or resi 85 and 2y9x or resi 94 and 2y9x or resi 200 and 2y9x or resi 256 and 2y9x or resi 259 and 2y9x or resi 263 and 2y9x or resi 295 and 2y9x
show licorice, cons_2y9x
color atomic, (cons_2y9x not elem C)
label n. CA and cons_2y9x, "%s-%s" % (resn, resi) 
select tropolone, 2y9x and resi 410
color atomic, tropolone
color yellow, (tropolone and elem C)

remove chain B or chain C or chain D or chain E or chain F or chain G or chain H or chain I or chain J

cealign 2y9x, 4j3p

hide (resi 2011 and 4j3p)
hide (resi 400 and 2y9x)
hide (resi 401 and 2y9x)

set sphere_scale, 0.5

set_view (\
     0.436669469,    0.420733839,   -0.795174658,\
     0.532359600,   -0.833371341,   -0.148599058,\
    -0.725194871,   -0.358431458,   -0.587891757,\
     0.000004590,   -0.000039697,  -80.084159851,\
    -7.154716492,  -25.432687759,  -32.313804626,\
    39.678813934,  120.495445251,  -20.000000000 )

@src/structural-visualizations/nicify.pml
set cartoon_transparency, 0.8
enable
#ray
#png manuscript/figures/structural-alignment-labels.png

hide labels
ray
png manuscript/figures/structural-alignment.png
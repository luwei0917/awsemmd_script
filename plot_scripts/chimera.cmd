open A.pdb
open B.pdb
morph start #0 frames 60
morph interpolate #1 frames 60
morph interpolate #1 frames 60
morph interpolate #0 frames 60
morph interpolate #0 frames 60
~ribbon #0
transparency 80 #1
morph movie
rainbow #2

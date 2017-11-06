

bg_color white
load final.pdb, final
# Hide all the atoms, then make the cartoon and display just the ribbons:
hide everything,final
cartoon automatic, final
show cartoon, final
set cartoon_fancy_helices=1

color red,ss h
color yellow ,ss s
color blue,ss l+''

spectrum count, rainbow_rev

ending
movie.add_roll(8.0,axis='x',start=1)
ending
movie.add_roll(8.0,axis='y',start=241)
movie.produce final.mpg,quality=90,quiet=0
ending

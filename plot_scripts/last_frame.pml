

bg_color white
load last.pdb, final
# Hide all the atoms, then make the cartoon and display just the ribbons:
hide everything,final
cartoon automatic, final
show cartoon, final
set cartoon_fancy_helices=1

color red,ss h
color yellow ,ss s
color blue,ss l+''

spectrum count, rainbow_rev

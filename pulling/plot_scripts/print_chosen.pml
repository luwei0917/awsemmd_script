bg_color white
load chosen.pdb, chosen
# Hide all the atoms, then make the cartoon and display just the ribbons:
hide everything,chosen
cartoon automatic, chosen
show cartoon, chosen
set cartoon_fancy_helices=1

color red,ss h
color yellow ,ss s
color blue,ss l+''

spectrum count, rainbow_rev

ray 2000,2000
png chosen.png

quit

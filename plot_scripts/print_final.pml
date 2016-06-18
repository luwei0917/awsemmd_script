bg_color white
load final.pdb, final
# Hide all the atoms, then make the cartoon and display just the ribbons:
hide everything, final
cartoon auto, final
show cartoon, final
set cartoon_fancy_helices=1

color red, ss h
color yellow, ss s
color green, ss l+''

ray 2000,2000
png final.png
quit

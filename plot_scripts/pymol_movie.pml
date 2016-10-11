load movie.pdb, mov
python

for x in cmd.get_names():

  cmd.create("mov", x, 1, -1)

python end
cmd.disable('all')
cmd.enable('mov',1)
cmd.hide("lines"     ,"all")
cmd.show("cartoon"   ,"mov")
load 2xov.pdb
cmd.hide("lines"     ,"2xov")
remove resn HOH
cmd.spectrum("count",selection="(mov)&*/CA")
sele resn DUM
cmd.show("spheres"   ,"sele")
mset 1 -400

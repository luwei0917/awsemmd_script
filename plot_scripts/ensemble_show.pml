bg_color white


load 1.pdb,base
python
for i in range(2,20):
  cmd.load(str(i)+".pdb")
  cmd.cealign("base",str(i))
  cmd.set("cartoon_transparency",0.5,str(i))
python end
orient animate=-1

# Hide all the atoms, then make the cartoon and display just the ribbons:
hide everything,all
cartoon automatic, all
show cartoon, all
set cartoon_fancy_helices=1
cmd.hide("lines"     ,"all")
util.color_deep("gray60", 'all')
util.color_deep("tv_red", 'base')

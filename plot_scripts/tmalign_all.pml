load result_all_atm.pdb
select a, c. A
cmd.show("cartoon"   ,"all")
util.color_deep("gray90", 'all')
util.color_deep("red", 'a')
cmd.hide("lines"     ,"all")
cmd.disable('a')
select native, chain B
select chosenframe, chain A

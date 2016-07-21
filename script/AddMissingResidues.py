from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# sample alignment.alignment
# >P1;1qjp
# structureX:1qjp:   1 :A:171 :A:::-1.00:-1.00
# APKDNTWYTGAKLGWSQ-------------HENKLGAGAFGGYQVNPYVGFEMGYDWLGRMPY-------AYKAQGVQLT
# AKLGYPITDDLDIYTRLGGMVWRADTYSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTN--------------G
# MLSLGVSYRFG*
#
# >P1;1qjp_fill
# sequence:::::::::
# APKDNTWYTGAKLGWSQYHDTGFINNNGPTHENKLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLT
# AKLGYPITDDLDIYTRLGGMVWRADTYSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDNG
# MLSLGVSYRFG*

class MyModel(automodel):
    def select_atoms(self):
        return selection(self.residue_range('18', '30'),
                         self.residue_range('64', '70'),
                         self.residue_range('146', '159'))

a = MyModel(env, alnfile = 'alignment.ali',
            knowns = '1qjp', sequence = '1qjp_fill')
a.starting_model= 1
a.ending_model  = 1

a.make()

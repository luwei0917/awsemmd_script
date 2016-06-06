import subprocess
import os
my_env = os.environ.copy()
my_env["PATH"] = "/home/wl45/python/bin:" + my_env["PATH"]
print(my_env)
subprocess.Popen("vmd", env=my_env)
os.system("vmd")
subprocess.run(["ls", "-l"])
# subprocess.run(["vmd", "-l"], shell=True)

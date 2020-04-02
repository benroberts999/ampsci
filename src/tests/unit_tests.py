#!/usr/bin/python3
import os
test_dir = "src/tests/"
filename = test_dir+"unit_tests.txt"
outfname = test_dir+"results.txt"
split = "#TEST"
output_seperator = "###########################################"
ifile = open(filename, "r")
#Read int the file, and strip any blank lines
content = os.linesep.join([s for s in ifile.read().splitlines() if s])

# If results file already exists, rename it to results_old (over-write)
if os.path.exists(outfname):
    old_name = outfname+"_old.txt"
    if os.path.exists(old_name):
        os.remove(old_name)
    os.rename(outfname,old_name)

import subprocess
import datetime
output = open(outfname,"a")
git_long = subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()
git_short = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()
branch = subprocess.check_output(['git', 'name-rev', '--name-only', 'HEAD']).decode('ascii').strip()
output.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")+"\n")
output.write("git branch: "+branch+"\n")
output.write("git commit: "+git_long+" ["+git_short+"]\n"+output_seperator+"\n\n")
output.close()

input_list = content.split(split)
#first item in list is blank (before first '#TEST#')
input_list.pop(0)
tmp_fname_input = "UNIT_TESTS_temp.DELETEME"
for input_item in input_list:
  input_file = open(tmp_fname_input,"w")
  output = open(outfname,"a")
  output.write(input_item+output_seperator+"\n")
  output.close()
  input_file.write(input_item)
  input_file.close()
  os.system("./diracSCAS "+tmp_fname_input+" |tee -a "+outfname)
  os.remove(tmp_fname_input)
  output = open(outfname,"a")
  output.write("\n"+output_seperator+"\n"+output_seperator+"\n")
  output.close()

import os
f = open('unrnames.txt','r')
f1 = open('unrmet.txt','w')
f1.write('[')
flines = f.readlines()
for line in flines[:-1]:
  rname = line.strip()
  cmd = "curl 'http://bigg.ucsd.edu/api/v2/universal/reactions/" + "%s'"%rname
  tmp = os.popen(cmd).read()
  f1.write(tmp + ',')
line = flines[-1]
rname = line.strip()
cmd = "curl 'http://bigg.ucsd.edu/api/v2/universal/reactions/" + "%s'"%rname
tmp = os.popen(cmd).read()
f1.write(tmp)
f1.write(']')
f1.close()
f.close()


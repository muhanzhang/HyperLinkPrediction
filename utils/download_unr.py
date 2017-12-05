import os
command = "curl 'http://bigg.ucsd.edu/api/v2/universal/reactions'"
unr = os.popen(command).read()
file = open('unr.txt', 'w')
file.write(unr)
file.close()            #to download all reactions from BIGG and save to rxns




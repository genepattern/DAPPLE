import sys
import os
import time
from datetime import date, timedelta

yesterday = date.today() - timedelta(1)
dateYesterday=yesterday.strftime('%m/%d/%Y')
if dateYesterday in open("/home/unix/rossin/dappleTracking/userTrackerEmailLog",'r').read():
	sys.exit()

#SENDMAIL = "/usr/sbin/sendmail" # sendmail location
#p = os.popen("%s -t" % SENDMAIL, "w")
file = open("/web/wwwprod/htdocs/mpg/dapple/emailTrackerV2",'r')

keep=[]
locations=[]
for line in file:
	x = line.split()[0:3]
	x = x + line.split("\t")[2:]
	if len(x)<6:
		continue
	if not x[0]=="liz.rossin@gmail.com" and not x[5]=="ERROR":
		if x[1]==dateYesterday:
			locations.append(line.split("\t")[8])
			keep.append(line)

out = open("/home/unix/rossin/dappleTracking/emailToSend",'w')
out.write("DAPPLE tracker from "+dateYesterday+"\n")
out.write(str(len(keep))+" jobs from "+", ".join(list(set(locations)))+"\n\n")
for item in keep:
	out.write(item.strip("\n")+"\n")

out.close()
#sts = p.close()

out = open("/home/unix/rossin/dappleTracking/userTrackerEmailLog",'w')
out.write(str(dateYesterday))
out.close()

os.system("chmod 777 /home/unix/rossin/dappleTracking/*")

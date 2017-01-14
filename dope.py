run = True
stupid = True
dumb = []
connum = 1
subnum = 97
while(run == True):
	var = raw_input("Contention = a, Subpoint = s, Evidence = d, Misc = f, Stop = end ")
		
	if(var is "a"):
		con = raw_input("What's the contention? ")
		dumb.append("Contention " + str(connum) + ": " + con)
		connum += 1
		stupid = False
		subnum = 97
		
		
	elif(var is "s" and stupid == False):
		sub = raw_input("What's the subpoint? ")
		dumb.append("Subpoint " + chr(subnum) + ": " + sub)
		subnum += 1
		
		
	elif(var is "d" and stupid == False):
		source = raw_input("What's the source? ")
		evi = raw_input("Summarize the card ")
		dumb.append("Evidence " + "(" + source + ")" + ": " + evi)

	elif(var is "f"):
		misc = raw_input("What do you want?")
		dumb.append("Comments: " + misc)
		
	elif(var == "end"):
		run = False
			
print dumb
for e in dumb:
	print(e)

from openpyxl import Workbook

wb = Workbook()

dest_filename = 'debate_flows1.xlsx'

ws1 = wb.active
ws1.title = "Debate Flow"

rick = 1
for e in dumb:
     d = ws1.cell(row = rick, column = 1, value = e)
     rick+=2

wb.save(filename = dest_filename)
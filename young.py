stupid = True
dumb = []
while(stupid == True):
	var = raw_input("Put stuff in: ")
	if(var == "end"):
		stupid = False
		break
	if((var[0]) is "/" and var[1].isdigit()):
		dumb.append("Contention " + var[1] + ": " + var[2:])
		break
	if((var[0]) is "."):
		dumb.append("Subpoint " + var[1] + ": " + var[2:])
		break
	else:
		dumb.append(var)
print dumb
for e in dumb:
	print(e)

##from openpyxl import Workbook

##wb = Workbook()

##dest_filename = 'debate_flows.xlsx'

##ws1 = wb.active
##ws1.title = "Debate Flow"

##rick = 1
##for e in dumb:
##     d = ws1.cell(row = rick, column = 1, value = e)
##     rick+=2

##wb.save(filename = dest_filename)
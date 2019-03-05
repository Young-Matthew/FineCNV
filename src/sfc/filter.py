# 

with open('/haplox/users/wangzy/04.FineCNV/09.SFClip/out/24_Combined.gDNA.dedup.sfc_breakpoints.location.txt') as file:
	print(file.readline(), end = '')
	for line in file:
		sep =  line.rstrip().split('\t')
		qualified = sum([True if i >10 else False for i in map(int, sep[2:])])
		if qualified > 0:
			print(line, end ='')



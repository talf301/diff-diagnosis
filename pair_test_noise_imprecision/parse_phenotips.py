import sys
import os

def get_codes(f):
	lines = f.readlines()
	codes = []
	for line in lines:
		if line.startswith('<span class="id"'):
			code = line[24:30]
			codes.append(code)
	return codes

def format_codes(codes, f):
	for code in codes:
		line = '9.9999\t9.9999\tOMIM:{}\tDISEASE NAME GOES HERE\tXXXX\t99999\n'.format(code)
		f.write(line)

if __name__ == '__main__':
	dirin = sys.argv[1]
	dirout = sys.argv[2]
	d = os.listdir(dirin)
	for fname in d:
		f = open(dirin + '/' + fname)
		fout = open(dirout + '/' + fname, 'w')
		omim_codes = get_codes(f)
		format_codes(omim_codes, fout)
		

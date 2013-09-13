
from sys import argv
from os import popen


f = popen('ls %s'%argv[1])
fns = f.read().split()
f.close()

for fn in fns:
    t = map(lambda(l): l.strip('\n'), open(fn).readlines())
    header = filter(lambda(l): l[0] == '#', t)
    txt = map(lambda(l): l.split('\t'), filter(lambda(l): l[0] != '#', t))
    
    newfn = fn.replace('.vcf', '.new.vcf')
    f = open(newfn, 'w')
    for l in header:
        f.write('%s\n'%l)
    for l in txt:
        if ':SS' in l[8]:
            f.write('%s\n'%('\t'.join(l)))
        else:
            newline = l
            newline[8] = newline[8] + ':SS'
            for i in range(9, len(newline)):
                newline[i] = newline[i] + ':.'
            f.write('%s\n'%('\t'.join(newline)))
    f.close()
    print 'finish %s'%fn
    

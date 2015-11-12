#!/usr/bin/env python

import sys

if len(sys.argv)!=3:
    print 'Usage: '+sys.argv[0]+' file.yoda directory'
    sys.exit(1)

begin = '# BEGIN'
end   = '# END'

hist = False
skip = 0

with open(sys.argv[1].rsplit('/',1)[-1].rsplit('.',1)[0]+'_'+sys.argv[2]+'.yoda','w') as fout:
    with open(sys.argv[1]) as fin:
        for line in fin:
            if skip > 0:
                skip -= 1
                continue
            elif hist:
                fout.write(line)
                if line[:len(end)]==end:
                    fout.write('\n')
                    hist = False
            elif line[:len(begin)]==begin:
                if '/'+sys.argv[2]+'/' in line:
                    name = line.rsplit('/',1)[1]
                    fout.write(
                        line.rsplit(' ',1)[0] +
                        ' /MC_HJETS_LH15/' + name +
                        'Path=/MC_HJETS_LH15/' + name
                    )
                    hist = True
                    skip = 1
    print fout.name


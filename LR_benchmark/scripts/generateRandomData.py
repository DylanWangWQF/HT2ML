#!/usr/bin/python

import random, sys, math

nFiles = 1
if len(sys.argv) < 4:
    print("usage: python generateRandomData.py filename d N [nFiles]")
    sys.exit(1)
else:
    filename = sys.argv[1]
    d = int(sys.argv[2])
    N = int(sys.argv[3])
    if len(sys.argv) > 4:
        nFiles = int(sys.argv[4])

MIN = 0
MAX = 2

valuesPerFile = int(math.ceil(float(N) / nFiles))

coeff = [random.uniform(0, 2) for i in range(d)]
for n in range(nFiles):
    if nFiles > 1:
        name = '%s_%d.dat' % (filename, n)
    else:
        name = '%s.dat' % filename

    f = open(name, 'w')

    if nFiles == 1 or n < nFiles - 1 or N % valuesPerFile == 0:
        nValues = valuesPerFile
    else:
        nValues = N % valuesPerFile

    f.write('%d %d\n' % (d, nValues))

    for i in range(nValues):
        val = [random.randint(MIN, MAX) for i in range(d)]
        # label = sum([coeff[i] * val[i] for i in range(d)])
        # label += random.gauss(0, 2);
        label = random.randint(MIN, MAX);

        for j in range(d):
            f.write('%d ' % val[j])
        f.write('%d\n' % label)
    f.close()


import sys

canonicalDict = {}
with open('ensCanonicalTranscript.txt') as f:
    for i in f:
        canonicalDict[i.split('\t')[0]] = i.split('\t')[1].strip()

# infile was annotated by annovar
infile = sys.argv[1]
of = open(sys.argv[2], 'w')
with open(infile) as f:
    header = f.readline()
    header = header.rstrip() + '\tcanonical_transcript\tpc\n'
    of.write(header)
    for line in f:
        info = line.strip().split('\t')
        if info[8] in ['nonsynonymous SNV']:
            if len(info[6].split(';')) == 1:
                aaList = info[9].split(';')
                for aa in aaList:
                    if aa.split(':')[1] == canonicalDict[aa.split(':')[0]]:
                        canonical = aa.split(':')[1]
                        pc = aa.split(':')[-1]
                        of.write('\t'.join(info) + '\t' + canonical + '\t' + pc + '\n')
                        break
                else:
                    canonical = aaList[0].split(':')[1]
                    pc = aaList[0].split(':')[-1]
                    of.write('\t'.join(info) + '\t' + canonical + '\t' + pc + '\n')
            else:
                for gene in info[6].split(';'):
                    tmpList = []
                    aaList = info[9].split(';')
                    for aa in aaList:
                        if aa.startswith(gene):
                            tmpList.append(aa)
                    if tmpList:
                        for aa in tmpList:
                            if aa.split(':')[1] == canonicalDict[aa.split(':')[0]]:
                                canonical = aa.split(':')[1]
                                pc = aa.split(':')[-1]
                                of.write('\t'.join(info) + '\t' + canonical + '\t' + pc + '\n')
                                break
                        else:
                            canonical = tmpList[0].split(':')[1]
                            pc = tmpList[0].split(':')[-1]
                            of.write('\t'.join(info) + '\t' + canonical + '\t' + pc + '\n')
of.close()

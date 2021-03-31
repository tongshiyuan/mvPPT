cds = {}
enstList = []
geneList = []
seqList = []
# Homo_sapiens.GRCh37.75.cds.all.fa can be download from
# ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz
with open('Homo_sapiens.GRCh37.75.cds.all.fa') as f:
    line = f.readline()
    seq = ''
    enst = line.lstrip('>').split()[0]
    gene = line.split()[3].lstrip('gene:')
    while True:
        line = f.readline()
        if line.startswith('>'):
            cds[gene] = cds.get(gene)
            if not cds[gene]:
                cds[gene] = {'enst': [], 'length': []}
            cds[gene]['enst'].append(enst)
            cds[gene]['length'].append(len(seq))
            enst = line.lstrip('>').split()[0]
            gene = line.split()[3].lstrip('gene:')
            seq = ''
        elif not line:
            break
        else:
            seq = seq + line.rstrip()
of = open('ensCanonicalTranscript.txt', 'w')
for k, v in cds.items():
    idx = v['length'].index(max(v['length']))
    of.write(k + '\t' + v['enst'][idx] + '\t' + str(max(v['length'])) + '\n')
of.close()

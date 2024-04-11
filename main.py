# !/usr/bin/env python3
import argparse
from Bio import AlignIO as AIO
from Bio import ExPASy
from Bio.Align import MultipleSeqAlignment as MSA
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from parsel import Selector
import re
import requests
import sys

re_baseid = re.compile(r'([^.]+)(\..+$)?')

class StrippedSequence:
    seq = ""
    stripped_region = -1

def parse_args(): 
    parser = argparse.ArgumentParser(prog=sys.argv[0], description="Convert stockholm MSA file ids to taxonomy ids.", add_help=True)
    parser.add_argument('-i', '--input', action='store', required=True, help="A PFAM stockholm file with canonical Uniprot identifiers for each sequence.")
    parser.add_argument('-o', '--output', action='store', default="out", help="Filename to write out a stockholm MSA file and unaligned FASTA file with NCBI taxonomy ids for each sequence.")
    parsed = parser.parse_args(sys.argv[1:])
    return(parsed)

if __name__ == '__main__':
    parsed=parse_args()
    file=parsed.input
    outfile_stockholm=parsed.output+'.sto'
    outfile_fasta=parsed.output+'.fasta'
    alignment=AIO.read(file, "stockholm")
    newalignment = []

    print('Translating record IDs...\n')
    ids = set()
    for record in alignment:
        accession = re_baseid.match(record.annotations['accession'])
        if (accession):
            try:
                uniprot_data = requests.get('https://www.uniprot.org/uniprot/{}.xml'.format(accession.group(1))).text
                sel = Selector(text=uniprot_data).css('dbReference[type="NCBI Taxonomy"]')
                sel_attrib = sel.attrib
                if 'id' in sel_attrib:
                    ncbi_taxid = sel_attrib['id']

                    # make NCBI taxid unique by concatenating counter
                    ncbi_taxid_revised = ncbi_taxid
                    if ncbi_taxid in ids:
                        counter = 1
                        ncbi_taxid_revised = ncbi_taxid + '_' + str(counter)
                        while ncbi_taxid_revised in ids:
                            counter += 1
                            ncbi_taxid_revised = ncbi_taxid + '_' + str(counter)
                    record.id = ncbi_taxid_revised
                    ids.add(ncbi_taxid_revised)
                    newalignment.append(record)
            except Exception as e:
                print(e)
                pass

    AIO.write(MSA(newalignment), outfile_stockholm,'stockholm')

    # write out as unaligned fasta
    x=list(SeqIO.parse(outfile_stockholm,'stockholm'))
    o = open(outfile_fasta, 'w')
    for i in x:
        o.write('>' + i.id[:i.id.find('/')] + '\n')
        s = str(i.seq)
        s = s.replace('-', '')
        o.write(s + '\n')

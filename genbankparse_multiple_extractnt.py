#!/usr/local/python/2.7.3/bin/python
from Bio import SeqIO
gb_file='21819_90A1.gbk'
nt=open ('90A1_nt.txt','a+')
import csv
with open ('90A1.csv','w') as aave:
    names=['locus_tag','note','gene','location','old_locus_tag','inference',\
           'product','nucleotide','translation','protein_id']
    writer=csv.DictWriter(aave,fieldnames=names,lineterminator='\n')
    writer.writeheader()
    for gb_record in SeqIO.parse(open(gb_file,'r'),'genbank'):
        for features in gb_record.features:
            qualifireslist=list(features.qualifiers)
            #print(qualifireslist)
            if features.type =='CDS':
                nucleotide=features.location.extract(gb_record).seq
                #print(type(nucleotide))
                locus_tag=features.qualifiers['locus_tag']
                if 'old_locus_tag' in qualifireslist:
                    old_locus_tag=features.qualifiers['old_locus_tag']
                else:
                    old_locus_tag=['no old locus tag']
                if 'note' in qualifireslist:
                    note=features.qualifiers['note']
                else:
                    note=['no note']
                if 'gene' in qualifireslist:
                    gene=features.qualifiers['gene']
                else:
                    gene=['no gene']
                if 'inference' in qualifireslist:
                    inference=features.qualifiers['inference']
                else:
                    inference=['no inference']
                if 'product' in qualifireslist:
                    product=features.qualifiers['product']
                else:
                    product=['no product']
                if 'protein_id' in qualifireslist:
                    protein_id=features.qualifiers['protein_id']
                else:
                    protein_id=['no protein id']
                if 'translation' in qualifireslist:
                    translation=features.qualifiers['translation']
                else:
                    translation=['no translation']
                

                locus_tagstring=','.join(locus_tag)
                oltstring=','.join(old_locus_tag)
                infstr=','.join(inference)
                productstr=','.join(product)
                protein_idstr=','.join(protein_id)
                transstr=','.join(translation)
                notestr=','.join(note)
                genestr=','.join(gene)
                neocleotidestr=''.join(nucleotide)
             
                writer.writerow({'locus_tag':locus_tagstring,'note':notestr,\
                                 'gene':genestr,'location':features.location,\
                                 'old_locus_tag':oltstring,'inference':infstr,\
                                 'product':productstr,'nucleotide':nucleotide,\
                                 'translation':transstr,'protein_id':protein_idstr})
                nt.write('>')
                nt.write(locus_tagstring)
                nt.write('\n')
                nt.write(neocleotidestr)
                nt.write('\n')
nt.close()
                


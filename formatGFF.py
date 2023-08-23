#!/bin/python
# coding: utf-8
## with python version > 3.4.5


import re, argparse, sys, os.path, pprint
from collections import defaultdict

##definition du nom de la class (formatGFF) a laquelle appartient l'objet
class formatGFF (object):
	def __init__(self):
		"""
		Global attributes of the object
		definition des attributs de l'objet
		"""
		self.parameters=''
		self.curated=''
		self.inputlines=defaultdict(lambda:defaultdict(list))
		self.input_ID=set()
		
	##definition des methodes appliquees a l'objet
	##main est la methode/fonction maitre
	def main(self):
		self.parseOptions()
		self.checkInputs()
		self.selectLines()
		self.parseGFF()
		self.updateGFF()


	def updateGFF(self):
		if self.parameters.ref==None:
			file_to_update=self.parameters.input
		else:
			file_to_update=self.parameters.ref

		output_file=open(re.sub(".gff3", "_updated.gff3", file_to_update), 'w')
		output_file.write("##gff-version\t3\n")

		with open(file_to_update, 'r') as f2:
			ref_ID_seen=set()
			for rline in f2.readlines():
				if rline.startswith('#')==False:
					li=rline.rstrip().split('\t')
					attributes=li[8].split(';')
					ref_attributes=self.getAttributes(attributes)
					#print(ref_attributes)
					
					##to deal with CDS created and join to existing CDS in artemis: 
					if 'Traes' in ref_attributes['ID']:
						refID=ref_attributes['ID'].split('.')[0]
					elif 'CDS:' in ref_attributes['ID']:
						refID=ref_attributes['Parent'].split('.')[0]
					else:
						sys.stderr.write("feature ID %s unknown\n" % str(ref_attributes['ID']))

					if refID in self.input_ID and refID not in ref_ID_seen:
						curated_lines_to_write=self.write_orderedGFF(inlines=self.inputlines[refID])
						for line in curated_lines_to_write:
							output_file.write(line+"\n")

						ref_ID_seen.add(refID)

					if refID not in self.input_ID:
						output_file.write(rline)

			##for new_genes created in artemis, i.e. genes in difference between the 2 sets "input_ID" and "ref_ID_seen":
			for new_gene in self.input_ID.difference(ref_ID_seen):
				curated_lines_to_write=self.write_orderedGFF(inlines=self.inputlines[new_gene])
				for line in curated_lines_to_write:
					output_file.write(line+"\n")


	def write_orderedGFF(self, inlines):
		ordered_lines=[]
		for f in ['gene', 'mRNA', 'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR']:
			if f in ['exon', 'CDS','five_prime_UTR', 'three_prime_UTR'] and f in inlines.keys():
				for sub_feature in inlines[f]:
					ordered_lines.append(sub_feature)
			elif f in ['gene', 'mRNA']:
				ordered_lines.append(inlines[f][0])

		return ordered_lines


	def parseGFF(self):
		##read gff that needs to be reordered and integrated to the original gff file
		####/!\ SCRIPT DOIT POUVOIR FONCTIONNER SI "manually curated" PRESENT POUR FEATURE GENE UNIQUEMENT
		with open(self.parameters.input, 'r') as f1:
			IDfeat_count=defaultdict(lambda:defaultdict(int))
			for iline in f1.readlines():
				if iline.startswith('#')==False:
					iline=iline.rstrip()
					li=iline.split('\t')
					attributes=li[8].split(';')
					
					input_attributes=self.getAttributes(attributes)
					##to deal with CDS created and join to existing CDS in artemis: 
					if 'Traes' in input_attributes['ID']:
						geneID=input_attributes['ID'].split('.')[0]
					elif 'CDS:' in input_attributes['ID']:
						geneID=input_attributes['Parent'].split('.')[0]
					else:
						sys.stderr.write(" feature ID unknown\n")

					feature=li[2]
					IDfeat_count[geneID][feature]+=1

					##select curated lines
					if any(s in input_attributes['Note'] for s in self.curated): 
						##format curated lines attributes and features IDs USING output from self.getAttributes
						list_attr=[]
						if feature in ['exon', 'CDS','five_prime_UTR', 'three_prime_UTR']:
							##format ID attribute
							if IDfeat_count[geneID]['mRNA']==0:
								list_attr.append('ID='+geneID+'.1.'+feature+str(IDfeat_count[geneID][feature]))
							else:
								list_attr.append('ID='+geneID+'.'+str(IDfeat_count[geneID]['mRNA'])+'.'+feature+str(IDfeat_count[geneID][feature]))
							##format other attributes
							for a in ['Parent','previous_id','Target','Note']:    ##supress 'Name' attribute that were previous ID given by MAGATT for exon/CDS features
								if a in input_attributes.keys():
									list_attr.append(a+'='+input_attributes[a])

						elif feature=='mRNA':
							list_attr.append('ID='+geneID+'.'+str(IDfeat_count[geneID]['mRNA']))
							for a in ['Name','Parent','previous_id','mapping','Note']:
								if a in input_attributes.keys():									
									list_attr.append(a+'='+input_attributes[a])

						elif feature=='gene':
							for a in ['ID','Name','previous_id','mapping','Note']:
								if a in input_attributes.keys():
									list_attr.append(a+'='+input_attributes[a])

						else:
							sys.stderr.write(" ERROR: unknown feature in line: %s \n" % str(iline))
						
						modifli='\t'.join(li[0:8])+'\t'+';'.join(list_attr)

						self.inputlines[geneID][feature].append(modifli)
						self.input_ID.add(geneID)

			##to have as many 'exon' than 'CDS', add exon features
			for g in self.inputlines:
				if len(self.inputlines[g]['CDS'])!=len(self.inputlines[g]['exon']):
					for cds in self.inputlines[g]['CDS']:
						exon=re.sub('CDS','exon', cds)
						self.inputlines[g]['exon'].append(exon)


	def getAttributes(self, attributes):
		dictio=defaultdict(str)
		output_attr=[]
		for a in attributes:
			if "Derives_from=" in a:
				a=re.sub('Derives_from=','ID=', a)
			(key,val)=a.split('=')
			##select one only attribute 'Note' in curated lines
			if key in ['Note','note']:
				if any(s in val for s in self.curated):
					dictio['Note']=val
			else:
				dictio[key]=val
		return dictio


	def selectLines(self):
		self.curated=self.parameters.select.rstrip('\n').split(',')


	def parseOptions(self):
		##use argparse to setup/get options from cmd line
		parser=argparse.ArgumentParser(description='Format a GFF3 after annotation curation and integrate changes into a reference GFF3 based on common geneIDs')
		parser.add_argument('-i', '--input', help='GFF3 file with curated annotation', required=True)
		parser.add_argument('-r', '--ref', help='reference GFF3 file to update')
		parser.add_argument('-s', '--select', help='comma separated strings to select curated lines from input GFF3 (default: manually,curated)', default='manually,curated')
		self.parameters=parser.parse_args()


	def checkInputs(self):
		"""
		check for input/output files
		"""
		sys.stdout.write('Check input files:\n')
		self.checkFile(file=self.parameters.input)
		if self.parameters.ref!=None:
			self.checkFile(file=self.parameters.ref)

	def checkFile (self,file):
		if os.path.isfile(file):
			sys.stdout.write(" file %s found\n" % str(file))
		else:
			sys.stderr.write(" ERROR: cannot find file %s \n" % str(file))
			sys.exit()


if __name__=='__main__':
		##creation d'un objet run de class formatGFF
		##execution de la methode/fonction main
		run=formatGFF()
		run.main()

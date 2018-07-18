#Aaron Till
#Genetic algorithm for producing a theoreticial evolutionary path between two strand of DNA

def main():
	origin_DNA = 'CGTCGTAATGACTGC'
	workingstrand = origin_DNA
	descendent_DNA = 'CGGCGAAACGATTGTTA'
	print('strand A', origin_DNA, '-> strand B', descendent_DNA)
	print()
	evolution(workingstrand, origin_DNA, descendent_DNA)



def Get_peptide(strand): # this is my code for converting a dna strand to its peptide sequence
	amino_acid = ''
	while len(strand) >= 3: 
		codon = strand[0:3]
		if codon == 'GCT' or codon == 'GCC' or codon == 'GCA' or codon == 'GCG':
			amino_acid += 'Ala'
			strand = strand[3:]
		if codon == 'CGT' or codon == 'CGC' or codon == 'CGA' or codon == 'CGG' or codon == 'AGA' or codon =='AGG': 
			amino_acid += 'Arg'
			strand = strand[3:]
		if codon == 'AAT' or codon =='AAC':
			amino_acid += 'Asn'
			strand = strand[3:]
		if codon == 'GAT' or codon =='GAC':
			amino_acid += 'Asp'
			strand = strand[3:]
		if codon == 'TGT' or codon == 'TGC':
			amino_acid += 'Cys'
			strand = strand[3:]
		if codon == 'CAA' or codon =='CAG':
			amino_acid += 'Gln'
			strand = strand[3:]
		if codon == 'GAA' or codon == 'GAG': 
			amino_acid += 'Glu'
			strand = strand[3:]
		if codon == 'GGT' or codon == 'GGC' or codon == 'GGA' or codon == 'GGG':
			amino_acid += 'Gly'
			strand = strand[3:]
		if codon == 'CAT' or codon == 'CAC':
			amino_acid += 'His'
			strand = strand[3:]
		if codon == 'ATT' or codon == 'ATC' or codon == 'ATA':
			amino_acid += 'Ile'
			strand = strand[3:]
		if codon == 'ATG':
			amino_acid += 'Met'
			strand = strand[3:]
		if codon == 'TTA' or codon =='TTG' or codon == 'CTT' or codon =='CTC' or codon == 'CTA' or codon == 'CTG':
			amino_acid += 'Leu'
			strand = strand[3:]
		if codon == 'AAA' or codon =='AAG':
			amino_acid += 'Lys'
			strand = strand[3:]
		if codon == 'TTT' or codon == 'TTC': 
			amino_acid += 'Phe'
			strand = strand[3:]
		if codon == 'CCT' or codon =='CCC' or codon =='CCA' or codon == 'CCG':
			amino_acid += 'Pro'
			strand = strand[3:]
		if codon == 'TCT' or codon == 'TCC' or codon == 'TCA' or codon == 'TCG' or codon == 'AGT' or codon == 'AGC':
			amino_acid += 'Ser'
			strand = strand[3:]
		if codon == 'ACT' or codon == 'ACC' or codon == 'ACA' or codon =='ACG': 
			amino_acid += 'Thr'
			strand = strand[3:]
		if codon == 'TGG':
			amino_acid += 'Trp'
			strand = strand[3:]
		if codon == 'TAT' or codon == 'TAC':
			amino_acid += 'Tyr'
			strand = strand[3:]
		if codon == 'GTT' or codon == 'GTC' or codon == 'GTA' or codon == 'GTG': 
			amino_acid += 'Val'
			strand = strand[3:]
		if codon == 'TAA' or codon == 'TGA' or codon == 'TAG':
			amino_acid += 'Stp'
			strand = strand[3:]
	return amino_acid


def evolution(workingstrand, origin_DNA, descendent_DNA): #this is my main program for imposing selection on the progenies
	history = {}
	history[0] = origin_DNA
	generation = 0
	peptide_change = {}
	peptide_change[0] = Get_peptide(origin_DNA)
	number_of_intermediates = 0
	temp_dictionary = {}
	minhamdistance = 200
	while workingstrand != descendent_DNA:
		generation += 1 
		progeny = replicate(workingstrand, generation) #calls the replicator function
		for i in range(len(progeny)):   #cheks list to see if any peptide matches descendent_DNA
			hammingtotal = hamming_distance(progeny[i], descendent_DNA) + (pep_dif(progeny[i], history[generation-1], workingstrand, history, generation))#two evaluation functions
			temp_dictionary[hammingtotal] = progeny[i] #if one does, then yay, we are done
			if hammingtotal < minhamdistance:
				minhamdistance = hammingtotal
		history[generation] = temp_dictionary[minhamdistance]
		peptide_change[generation] = Get_peptide(temp_dictionary[minhamdistance])
		if Get_peptide(temp_dictionary[minhamdistance]) != Get_peptide(origin_DNA) and  Get_peptide(temp_dictionary[minhamdistance]) != Get_peptide(descendent_DNA): 
			number_of_intermediates += 1
		workingstrand = temp_dictionary[minhamdistance]
	print()
	print('This dictionary shows which strand was selected within its generation')
	print(history)
	print()
	print('This dictionary shows the change in peptide sequence, there were', number_of_intermediates, 'intermediate peptides')
	print(peptide_change)
	return



def replicate(workingstrand, generation):
	progeny = [workingstrand]  #creates a list of copies itself with every possible mutation at each point. 
	for i in range(len(workingstrand)):
		temporarystrand = workingstrand[0:i] + 'C' + workingstrand[i+1:len(workingstrand)] #SUBSTITUITIONS here and below
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + 'A' + workingstrand[i+1:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + 'G' + workingstrand[i+1:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + 'T' + workingstrand[i+1:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand) 
		temporarystrand = workingstrand[0:i] + 'C' + workingstrand[i:len(workingstrand)] # INSERTIONS here and below
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + 'A' + workingstrand[i:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + 'G' + workingstrand[i:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + 'T' + workingstrand[i:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] +  workingstrand[i+1:len(workingstrand)] #DELETIONS here and below
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + workingstrand[i+1:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] +  workingstrand[i+1:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand)
		temporarystrand = workingstrand[0:i] + workingstrand[i+1:len(workingstrand)]
		if temporarystrand != workingstrand:
			progeny.append(temporarystrand) 
	#print()		
	#print('Progeny of', generation, 'consists of', progeny)
	return progeny


def hamming_distance(string1, string2): #selection for reaching descendent DNA
	ham_distance = 0
	for i in range(min(len(string1),len(string2))):
		if string1[i] != string2[i]:
			ham_distance += 1
	ham_distance += 2*(max(len(string1),len(string2)) - min(len(string1),len(string2))) #so different length is penalized
	return ham_distance


def pep_dif(string1, string2, workingstrand, history, generation): #selection against forming intermediates
	tracker = 0
	distance = 0
	pep1 = Get_peptide(string1)
	pep2 = Get_peptide(string2)
	for i in range(len(string1)//3*3): #this is so the tails don't form their own mismatches
		if pep1[i*3-3:i*3] == 'Stp': #strong selection against premature stop codon
			distance += 10
		if pep1[i*3-3:i*3] != pep2[i*3-3:i*3]:
			for i in range(len(history)-1):
				if history[i] == workingstrand: #preventing same selection as well as loops
					tracker = 0 #if hamming must be forced, gets rid of selection against intermediates 
				else: 
					tracker = 1
			distance += tracker
	return distance





main()

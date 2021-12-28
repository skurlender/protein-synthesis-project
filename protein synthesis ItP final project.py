'''
Sam Kurlender
skurlender@wesleyan.edu
COMP 112, Section 1
2018.11.28
'''

import os
import os.path

def import_DNA(file_path):
    '''
    signature: string -> string
    this function takes a string input, and if the input exists as a file, returns a string comprised of all the letters in the file
    '''
    dna_string = ""
    try:
        f = open(file_path, 'r')
        with open(file_path) as f:
            dna_list = f.read().splitlines()
        for letter in dna_list:
            dna_string = dna_string + (str(letter.upper()))
        return dna_string
    except FileNotFoundError:
        print("Sorry, " + file_path + " does not exist")


def is_valid_DNA(isvalid_string):
    '''
    singature: string -> boolean
    this function checks to see if the string from the file is a valid DNA string, meaning that A, T, C, and G are the only letters in it
    '''
    i = 0
    while(i < len(isvalid_string)):
        if(isvalid_string[i] != 'A' and isvalid_string[i] != 'T' and isvalid_string[i] != 'C' and isvalid_string[i] != 'G' and isvalid_string[i] != 'a' and isvalid_string[i] != 't' and isvalid_string[i] != 'c' and isvalid_string[i] != 'g'):
            return False
        else:
            i = i + 1
    return True

def RNA_complement(dna_string):
    '''
    signature: string -> string
    this function takes the valid DNA string and returns the RNA complement of the DNA string
    '''
    newstring = ""
    i = 0
    while(i<(len(dna_string))):
        if(dna_string[i] == 'A' or dna_string[i] == 'a'):
            newstring = newstring + 'U'
        elif(dna_string[i] == 'T' or dna_string[i] == 't'):
            newstring = newstring + 'A'
        elif(dna_string[i] == 'C' or dna_string[i] == 'c'):
            newstring = newstring + 'G'
        elif(dna_string[i] == 'G' or dna_string[i] == 'g'):
            newstring = newstring + 'C'
        i = i + 1
    return newstring


def find_start_position(rna_string):
    '''
    signature: string -> int
    this function takes the RNA complement and returns the index right after the letters 'AUG' appear
    '''
    i = 0
    while(i<len(rna_string)-1):
        if(rna_string[i] == 'A' and rna_string[i+1] == 'U' and rna_string[i+2] == 'G'):
            return i+3
        else:
            i = i + 1

def build_protein(rnastrand, index):
    '''
    signature: string, int -> string
    this function takes the RNA complement and the start position and returns the protein sequence that the RNA complement is made of
    '''
    codon_table = {"UUU":"F","UUC":"F","UUA":"L","UUG":"L","CUU":"L","CUC":"L","CUA":"L","CUG":"L",
                   "AUU":"I","AUC":"I","AUA":"I","AUG":"M","GUU":"V","GUC":"V","GUA":"V","GUG":"V",
                   "UCU":"S","UCC":"S","UCG":"S","CCU":"P","CCC":"P","CCA":"P","CCG":"P","ACU":"T",
                   "ACC":"T","ACA":"T","ACG":"T","GCU":"A","GCC":"A","GCA":"A","GCG":"A","UAU":"Y",
                   "UAC":"Y","CAU":"H","CAC":"H","CAA":"Q","CAG":"Q","AAU":"N","AAC":"N","AAA":"K",
                   "AAG":"K","GAU":"D","GAC":"D","GAA":"E","GAG":"E","UGU":"C","UGC":"C","UGG":"W",
                   "CGU":"R","CGC":"R","CGA":"R","CGG":"R","AGU":"S","AGC":"S","AGA":"R","AGG":"R",
                   "GGU":"G","GGC":"G","GGA":"G","GGA":"G","GGG":"G"}
  
    newstring = ""
    i = index
    while (i<len(rnastrand)):
        if(rnastrand[i:i+3] in codon_table):
            newstring = newstring + codon_table[rnastrand[i:i+3]]
        if(rnastrand[i:i+3] == "UAA" or rnastrand[i:i+3] == "UAG" or rnastrand[i:i+3] == "UGA"):
            return newstring
# these three strings are stop codons. if they appear in the string, we stop iterating through the string and just return the new string
        i = i + 3
    return newstring


def DNA_to_protein(file_path):
    '''
    signature: string -> string
    this function ties together all the code from the previous parts into one digestible function
    '''
    dnastring = import_DNA(file_path)
    if(is_valid_DNA(dnastring) == True):
        RNAstring = RNA_complement(dnastring)
        ind = find_start_position(RNAstring)
        proteins = build_protein(RNAstring,ind)
        return proteins
    else:
        string = "This is not a valid string of DNA, try again"
        return string
    

def DNAstrings(bigDNA, smallDNA):
    '''
    signature: string, string -> tuple(boolean, int)
    this function sees if the smallDNA string repeats itself twice in the bigDNA string. if it does, it returns a tuple of a boolean, either true or false, and the index of the beginning of the repetition  
    '''
    i = 0
    while (i<(len(bigDNA)-1)):
        if((bigDNA[i:i+3] == smallDNA) and (bigDNA[i+3:i+6] == smallDNA)):
            return (True,i)
        i = i + 1
    return (False,i)

def DNAstringrep(bigDNA, smallDNA):
    '''
    signature: string,string -> int
    this function returns the number of times the smallDNA string occurs in the bigDNA string
    '''
    count = 0
    it = 0
    z = DNAstrings(bigDNA, smallDNA)
    r = z[1]
    if(z[0] == True):
        while(it<len(bigDNA)-1):
            if(bigDNA[r:r+3] == smallDNA):
                r = r + 3
                count = count + 1
            it = it + 1
    return count
        
def main():
    '''
    this is the main function. it gives a testable condition for all of the above code
    '''
    file_path = "proteintest.txt"
    x = DNA_to_protein(file_path)
    print("I took a DNA string, translated it to an RNA string, and then returned the amino acid sequence that makes up the string of RNA. The sequence is " + x)
    bigDNA = "AGGACATGCATGCCTAGGAGGAGGAGGAGGAGGTAG"
    smallDNA = "AGG"
    t = DNAstringrep(bigDNA, smallDNA)
    print("I took two strings of DNA, a large string and a small string, and returned the number of times the small DNA string appears in the large DNA string. It appears " + str(t) + " times")    
    



    

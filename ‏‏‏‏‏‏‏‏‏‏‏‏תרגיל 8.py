#Shira Cohen 211485198
#Motive finding - Gibbs sampler algorithm

import random

#This function Randomly chooses starting positions
#and form the set of k-mers associated with these starting positions
def FindMotifs(dna, k):
    motifs=[]
    for i in dna:
        start=random.randint(0,len(i)-k)
        motifs+=[i[start:start+k]]
    return motifs

#This function generates a position specific counting matrix
#without adding psuidocount.
def CountOccurances(motifs, k):
    counter=[]
    for i in range(k): #Initializes the matrix
        counter+=[[0,0,0,0]]
    for i in motifs:
        for j in range(k):
            if i[j]=='A':
                counter[j][0]+=1
            elif i[j]=='T':
                counter[j][1]+=1
            elif i[j]=='C':
                counter[j][2]+=1        
            elif i[j]=='G':
                counter[j][3]+=1
    return counter

#This function counts the number of mismatches between nucleotides in the list "motifs"
#and returns the sum of the mismatch counts as the score.
def CalcScore(motifs):
    score=0
    counter=CountOccurances(motifs,k)
    for i in counter:
        for j in i:
            score+=j
        score-=max(i)
    return score

#This function generates a matrix containing the probabilities of PSSM 
def FindProfile(motifs):
    probabilities=[]
    for i in range(k): #Initializes the matrix
        probabilities+=[[0,0,0,0]]
    counter=CountOccurances(motifs,k)
    x=counter[0][0]+counter[0][1]+counter[0][2]+counter[0][3]+4 #Number of sequences + pseudocounts
    for i in range(len(probabilities)):
        for j in range(len(probabilities[i])):
            probabilities[i][j]+=(float)(counter[i][j]+1)/x
    return probabilities

#This function computes the probabilities of all k-mers in the removed sequence        
def ProfileProb(sequence, profile):
    probabilities=[]
    for i in range(len(sequence)-k+1):
        kmer=sequence[i:i+k]
        x=1
        for j in range(len(kmer)):
            if kmer[j]=='A':
                x*=profile[j][0]
            elif kmer[j]=='T':
                x*=profile[j][1]
            elif kmer[j]=='C':
                x*=profile[j][2]
            elif kmer[j]=='G':
                x*=profile[j][3]
        probabilities+=[x]
    return probabilities

#This function implements the Gibbs sampling algorithm
def GibbsSampler(dna, k, N):
    motifs=FindMotifs(dna,k)
    BestMotifs = motifs
    for i in range(N):
        x=random.randint(0,len(dna)-1)
        sequence=dna[x]
        del motifs[x]
        profile=FindProfile(motifs)
        profileProb=ProfileProb(sequence, profile)
        maxProb=max(profileProb)
        index=profileProb.index(maxProb)
        motifs.insert(x,sequence[index:index+k])
        score=CalcScore(motifs)
        if score<CalcScore(BestMotifs):
            BestMotifs=motifs

    return BestMotifs,score

#This function gets a filename and returns a list with the DNA sequences.
def GetDNA(filename):
    seq_file=open(filename,'r')
    dna=[]
    for line in seq_file:
        dna+=[line.strip()]
    seq_file.close()
    return dna

#This function repeats the whole process in order to reach convergence
def repeatGibbsSampler(dna, k, N, repeats):
    bestMotifs = FindMotifs(dna,k)    
    bestScore = CalcScore(bestMotifs)
    
    for i in range(repeats):
        (motifs, score) = GibbsSampler(dna, k, N)
        if score < bestScore:
            bestMotifs = motifs
            bestScore = score
    return bestMotifs
    
dna=GetDNA('DNA_sequences.txt')
k=4
N=100
repeats=100
print "The motifs are:"
print repeatGibbsSampler(dna, k, N, repeats)

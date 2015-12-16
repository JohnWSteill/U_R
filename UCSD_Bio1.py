ref = 'ACGT'
comp = dict([('A','T'), ('T','A'), ('C','G'), ('G','C')])
bases = ['U', 'C', 'A', 'G']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
amino_mass = {'A':   71.03711, 'C':   103.00919, 'D':   115.02694,
              'E':   129.04259,'F':   147.06841, 'G':   57.02146,
              'H':   137.05891,'I':   113.08406, 'K':   128.09496,
              'L':   113.08406,'M':   131.04049, 'N':   114.04293,
              'P':   97.05276, 'Q':   128.05858, 'R':   156.10111,
              'S':   87.03203, 'T':   101.04768, 'V':  99.06841,
              'W':  186.07931, 'Y':   163.06333} 
codon_table = dict(zip(codons, amino_acids))
codon_table_inv = {}
for key in set(codon_table.values()):
    codon_table_inv[key] = []
for key in codon_table.keys():
    codon_table_inv[codon_table[key]].append(key)
    
def ros_rmna(text):
    text += '*'
    count = 1
    for el in text:
        count *= len(codon_table_inv[el])
        count %= 1000000
    return count

from math import factorial as fac
def ros_lia(k,N):
    # Given: Two positive integers k N . In this problem, we begin with Tom,
    # who in the 0th generation has genotype Aa Bb. Tom has two children in
    # the 1st generation, each of whom has two children, and so on. Each
    # organism always mates with an organism having genotype Aa Bb.
    #
    # Return: The probability that at least N Aa Bb organisms will belong to
    # the k-th generation of Tom's family tree (don't count the Aa Bb mates at
    # each level). Assume that Mendel's second law holds for the factors.
    #   Aa  Aa  100
    # AA  Aa  aa  25-50-25 x Aa
    # AA Aa  | AA Aa aa |  Aa aa  25/2 25/2 | 25/2 50/2 25/2 | 25/2 25/2
    #  AA Aa aa : 25 50 25
    p = .25 # 
    pop = 2**k
    prob_atLeast = 0
    for n in range(N,pop+1):        
        prob_exact = 1.* fac(pop)/fac(n)/fac(pop-n) * p**n * (1-p)**(pop-n)
        prob_atLeast += prob_exact
        print (n, prob_exact, prob_atLeast)
    return prob_atLeast
    
def PatternCount(Text, Pattern):
    count = srchWindowBegin = 0        
    while True:
        srchWindowBegin = Text.find(Pattern, srchWindowBegin) + 1
        if srchWindowBegin > 0:
            count+=1
        else:
            return count

def countNucleotides(Text):
    counts = [ PatternCount(Text,char) for char in 'ACGT']
    print (' '.join([str(i) for i in counts]))
    return counts

def countNucleotidesFromFile(filename):
    with open(filename) as f:
        c = countNucleotides(f.readline().strip())
        
def DNAtoRNA(text):
    return text.replace('T','U')

def printDNAtoRNA(filename):
    with open(filename) as f:
        print (DNAtoRNA(f.readline().strip()))
        

def FibNK(n,k):
    if n == 0 : return 1
    if n == 1: return 1
    F = [ 1 for i in range(n)]
    for i in range(2,n):
        F[i] = F[i-1] + k * F[i-2]
    return F[-1]

def FibNM(n,m):
    if n == 1:
        return 1
    pop = [ [0 for i in range(m)]]
    pop[0][0] = 1
    for i in range(1,n):
        pop.append([0 for j in range(m)])
        pop[i][1:m] = pop[i-1][0:m-1]
        pop[i][0] = sum(pop[i-1][1:m])
    return sum(pop[-1])
    
    n = n+1 # "after nth month"
##    Return: The total number of pairs of rabbits that will remain after the n-th month if all
##    rabbits live for m months.
    if n == 0:
        return 1
    elif n == 1:
        return 1
    else:
        fibnm = [1 for i in range(n)]
        for i in range(2,n):
            fibnm[i] = fibnm[i-1] + fibnm[i-2]
            if i >= m:
                fibnm[i] -= fibnm[i-m] 
    return fibnm
                
                    

def printFib(filename):
    with open(filename) as f:
        [n,k] = [int(el) for el in f.readline().strip().split(' ')]
    print (FibNK(n,k))

def ros_iprb(k,m,n):
    # Given: Three positive integers k, m, and n, representing a population
    # containing k+m+n organisms: k individuals are homozygous dominant for a
    # factor, m are heterozygous, and n are homozygous recessive.

    # Return: The probability that two randomly selected mating organisms will
    # produce an individual possessing a dominant allele (and thus displaying
    # the dominant phenotype). Assume that any two organisms can mate.

    count = float(k+m+n);
    
    prob_kk = k/count * (k-1)/(count-1) if (k>1) else 0.
    prob_km = 2*k/count *  m/(count-1) if (k >0 and m>0) else 0.
    prob_kn = 2*k/count *  n/(count-1) if (k >0 and n>0) else 0.
    prob_mm = m/count * (m-1)/(count-1) if (m>1) else 0.
    prob_mn = 2*m/count *  n/(count-1) if (m >0 and n>0) else 0.
    prob_nn = 1 - sum([prob_kk,prob_km,prob_kn,prob_mm,prob_mn])

    return sum([prob_kk,prob_km,prob_kn]) + .75*prob_mm + .5*prob_mn
    

def printRI(filename):
    with open(filename) as f:
        [k,m,n] = [int(el) for el in f.readline().strip().split(' ')]
    print (ros_iprb(k,m,n))

def ros_gc(filename):
    (maxGcVal,maxkey) = (-1,"")
    
    def updateMM(k,t,maxV,maxK):
        tc = countNucleotides(text)
        gcVal= float(tc[1]+tc[2])/float(sum(tc))
        if gcVal > maxV:
            maxV = gcVal
            maxK = key
        return (maxV,maxK)
   
    f = open(filename)
    key = f.readline().strip()[1:]
    text = ''
    for line in f:
        nextline = line.strip()
        if (nextline[0]=='>'):
            (maxGcVal,maxkey) = updateMM(key,text,maxGcVal,maxkey)
            key = nextline[1:]
            text = ''
        else:
            text = text + nextline
    (maxGcVal,maxkey) = updateMM(key,text,maxGcVal,maxkey)        
    print (maxkey)
    print (maxGcVal*100)



def readFasta(fastaStream):
    key = fastaStream.__next__()[1:]
    text = ''
    for line in fastaStream:
        nextline = line.strip()
        if (nextline[0]=='>'):
            yield (key, text)
            key = nextline[1:]
            text = ''
        else:
            text = text + nextline
    yield (key, text)

def ros_mrpt(inp, pat = ['N','{P}','[ST]','{P}']):
    tests = []
    for req in pat:
        if len(req) == 1:
            tests.append(lambda x: x==pat)
        elif req[0] == '{':
            tests.append(lambda x: not any([x==pat[i]
                                            for i in range(1,len(req)-1)]))
        elif req[0]=='[':
            tests.append(lambda x: any([x==pat[i]
                                        for i in range(1,len(req)-1)]))
    for prot in inp:
        url = 'http://www.uniprot.org/uniprot/'+prot[:6]+'.fasta'
        
    
        

def ros_grph(filename, overlap_k = 3):
    sufx = {}
    prfx = {}
    adjL = []
    with open(filename) as f:
        for (key,text) in readFasta(f):
            sufx[key]= text[-overlap_k:]        
            prfx[key]= text[:overlap_k]
    for p in iter(sorted(prfx.items(), key=lambda x: x[1])):
        for s in iter(sorted(sufx.items(), key=lambda x: x[1])):
            if p[1]== s[1] and p[0] != s[0] :
                adjL.append((s[0],p[0]))
    for el in adjL:        
        print (' '.join([el2 for el2 in el]))
    return adjL
             
def getConsensusFromProfile(profile):
    cons = []
    for j_col in range(len(profile[0])):
        pvals = [profile[i_row][j_col] for i_row in range(4)]
        cons.append(ref[pvals.index(max(pvals))])
    return ''.join(cons)
    
def ros_cons(filename):
    dna = list(el[1] for el in readFasta(filename))
    profile = getProfileFromMotif(dna)
    return getConsensusFromProfile(profile)

    
def ros_prot(text):
    return ''.join([codon_table[text[i:i+3]] for i in range(0,len(text),3)])

def ros_prtm(prot_text):
    return sum(amino_mass[el] for el in prot_text)
    
def strToKmerInd(kmer,k):
    return sum([ref.find(kmer[i])*(4**(k-i-1)) for i in range(k)])

def kmerIndToStr(val,k):
    retStr = ['-' for i in range(k)]
    for i in range(k):
        [val, retStr[k-i-1]] = [val/4, ref[val%4]]
    return ''.join(retStr)

def FreqWords(text,k): # Using Index
    count = buildFreqArray(text,k)
    maxCount = max(count)
    inds = [ ind for ind in range(len(count)) if (count[ind] == maxCount)]
    FrequentPatterns = [kmerIndToStr(i) for i in inds]
    return FrequentPatterns
    
def revComp(Text):
    return ''.join([comp[z] for z in Text[::-1]])

def printRevComp(filename):
    with open(filename) as f:
        print(revComp(f.readline().strip()))
        
    
def findOcc(Text, Pattern):
    inds = []
    srchWindowBegin = 0        
    while True:
        srchWindowBegin = Text.find(Pattern, srchWindowBegin) + 1
        if srchWindowBegin > 0:
            inds.append(srchWindowBegin-1)
        else:
            return inds

def printFO(filename):
    with open(filename) as f:
        [text, pat] = [el.strip() for el in f.readlines()]
    inds = findOcc(text, pat)
    print(' '.join([str(i+1) for i in inds]))

def findInds(myList, val):
    inds = []
    srchWindowBegin = 0        
    while True:
        try:
            srchWindowBegin = myList.index(val, srchWindowBegin) + 1
        except:
            return inds
        if srchWindowBegin > 0:
            inds.append(srchWindowBegin-1)
        

def buildFreqArray(text,k):
    count = [ 0 for i in range(4**k)]
    for i in range(len(text)-k):
        count[strToKmerInd(text[i:i+k],k)]+=1
    return count

def runBFA(filename):
    with open(filename) as f:
        [text,k] = [l.strip() for l in f.readlines()];
        k = int(k)
    printCommaSepList(buildFreqArray(text,k))
    

def getClumpKmers(Genome,k,t,L):
    clumpInds = [False for i in range(4**k)]
    freqArray = buildFreqArray(Genome[0:L],k)
    for i in range(len(freqArray)):
        clumpInds[i] = (freqArray[i] >= t)
    for i in range(1,len(Genome)-L):
        freqArray[strToKmerInd(Genome[i-1:i-1+k],k)] -= 1
        endKmerInd = strToKmerInd(Genome[i+L-k:i+L],k)
        freqArray[endKmerInd] += 1
        if (freqArray[endKmerInd] >= t):
            clumpInds[endKmerInd] = True
    return [ kmerIndToStr(i,k) for i in range(4**k) if clumpInds[i]]
  
def readClumpInp(filename):
    with open(filename) as f:
        [Genome, kLt] = [l.strip() for l in f.readlines()]        
        [k,L,t] = [int(i) for i in kLt.split(' ')]
    return {'Genome':Genome, 'k':k, 'L':L, 't':t}

def getSkew(Genome):
    skew = [0 for i in range(len(Genome)+1)]
    for i in range(len(Genome)):
        skew[i+1] =  skew[i] - int(Genome[i]=='C') + int(Genome[i]=='G');
    return skew

def printSkew(filename):
    with open(filename) as f:
       printSpaceSepList(getSkew(f.readline().strip()))

def printSpaceSepList(myList):
    print(' '.join([str(el) for el in myList]))
    
def printMinSkew(filename):
    with open(filename) as f:
        skew = getSkew(f.readline().strip())
    printSpaceSepList(findInds(skew,min(skew)))     

def getHammingDist(str1,str2):
    return sum([1 for (char1,char2) in zip(str1,str2) if (char1 != char2)])

def printHammingDist(filename):
    with open(filename) as f:
        [g1,g2] = [l.strip() for l in f.readlines()]
    print(getHammingDist(g1,g2))

def getApproxMatchInds(g,pat,maxD):
    inds = [];
    patL = len(pat)
    for i in range(len(g)-len(pat)+1):
        if (getHammingDist(g[i:i+patL],pat) <= maxD):
            inds.append(i)
    #print ' '.join([str(i) for i in inds])
    return inds

def doGetApproxMatchInds(filename):
    with open(filename) as f:
        [pat,g,maxD] = [l.strip() for l in f.readlines()]
    return len(getApproxMatchInds(g,pat,int(maxD)))

from itertools import combinations, product 
def getNeighbors(g,d):
    yield(g)
    k = len(g)
    g = list(g)
    neighbors = [g]
    alpha = list('ACGT')
    for i_d in range(1,d+1):
        for indices in list(list(i) for i in combinations(range(k),i_d)):
            possMutations = [ [alpha[i] for i in range(4)
                               if (alpha[i] != g[ind])]
                              for ind in indices]
            for mutation in product(*possMutations):            
                gm = g[:];
                for ind in range(i_d):
                    gm[indices[ind]] = mutation[ind]
##                    print (indices, possMutations, mutations,
##                           gm, ind, indices[ind])
                #neighbors.append(gm)
                yield ''.join(gm)   
    #return [''.join(nb) for nb in neighbors]
        
        

def getFreqWordsWithMismatches(g,k,d):
    count = [ 0 for i in range(4**k)]
    for substr in [g[i:i+k] for i in range(len(g)-k+1)]:
        for pat in getNeighbors(substr,d):
            count[strToKmerInd(pat,k)] += 1
    maxcount = max(count)
    return [kmerIndToStr(i,k) for i in range(4**k) if (count[i]==maxcount)]
    
def printGFWWM(filename):
    with open(filename) as f:
        [g, kd] =[l.strip() for l in f.readlines()]
        [k,d]= [int(el) for el in kd.split(' ')]
    print(' '.join(getFreqWordsWithMismatches(g,k,d)))

def getFreqWordsWithMismatchesAndRevComp(g,k,d):
    count = [ 0 for i in range(4**k)]
    for substr in [g[i:i+k] for i in range(len(g)-k+1)]:
        for pat in getNeighbors(substr,d):
            count[strToKmerInd(pat,k)] += 1
            count[strToKmerInd(revComp(pat),k)] += 1
    maxcount = max(count)
    return [kmerIndToStr(i,k) for i in range(4**k) if (count[i]==maxcount)]
    
def printGFWWMARC(filename):
    with open(filename) as f:
        [g, kd] =[l.strip() for l in f.readlines()]
        [k,d]= [int(el) for el in kd.split(' ')]
    print(' '.join(getFreqWordsWithMismatchesAndRevComp(g,k,d)))

def getkmer(dna,k):
    for text in dna:
        for ikmer in [text[i:i+k] for i in range(len(text)-k+1)]:
              yield ikmer

def getkmer_p(dna,k,d):
    for kmer in getkmer(dna,k):
        for kmer_p in getNeighbors(kmer,d):
            yield kmer_p

def isPatInAll(pat_p,dna,d):
    k = len(pat_p)
    for text in dna:
        inText = False
        for kmer in [text[i:i+k] for i in range(len(text)-k+1)]:
            inText = (inText or (getHammingDist(pat_p,kmer)<=d))
            #print (pat_p, kmer, inText)
        if (not inText):
            return False
    return True
        
def getMotifsBruteForce(dna,k,d):
    motifPats = []
    notChecked = [True for i in xrange(4**k)]
    for pat_p in getkmer_p(dna,k,d):
        if (notChecked[strToKmerInd(pat_p,k)]):
            notChecked[strToKmerInd(pat_p,k)]=False
            if isPatInAll(pat_p,dna,d):
                yield pat_p 
    
def printGMBF(filename):
    with open(filename) as f:
        [k,d] = [int(el) for el in f.readline().strip().split(' ')]
        dna = [el.strip() for el in f.readlines()]
    #isPatInAll('ATT',dna,1)
    print(' '.join([el for el in getMotifsBruteForce(dna,k,d)]))

from math import log                 
def getEntropy(filename):
    with open(filename) as f:
        motifs = [el.strip().upper().split('   ') for el in f.readlines()]
    ncol = len(motifs[0])
    profile = [[0 for i in range(ncol)] for j in range(4)]
    for i_row in range(len(motifs)):
        for j_col in range(ncol):
            profile[ref.find(motifs[i_row][j_col])][j_col]+=1
    profEntropy = 0;
    for j_col in range(ncol):
        tot = float(sum(profile[i][j_col] for i in range(4)))
        for i in range(4):
            profile[i][j_col] /= tot
            if (profile[i][j_col] > 1e-6):
                profEntropy -= profile[i][j_col] * log(profile[i][j_col],2)
    return profEntropy

def getAllKmers(k):
    for i in xrange(4**k):
        yield kmerIndToStr(i,k)

def getSubKmers(text,k):
    for i in xrange(len(text)-k+1):
        yield text[i:i+k]

def getMotifDist(kmer,dna,k):
    d = 0
    for text in dna:
        d_i = 1e10
        for subKmer in getSubKmers(text,k):
            d_i = min(d_i,getHammingDist(kmer,subKmer))
        d += d_i
    return d
            
def getMedianString(dna,k):
    distance = 1e10
    for kmer in getAllKmers(k):
        motifDist = getMotifDist(kmer,dna,k)
        if motifDist < distance:
            distance = motifDist
            median = [kmer]
        elif motifDist == distance:
            median.append(kmer)
    return median

def printGMS(filename):
    with open(filename) as f:
        k = int(f.readline().strip())
        dna = [el.strip() for el in f.readlines()]
    print(getMedianString(dna,k))

def getProfProb(kmer,profile,k):
    prob = 1.
    for j_col in xrange(k):
        i_row = ref.find(kmer[j_col])
        prob *= profile[i_row][j_col]
    return prob

def getProfileMostProbKmer(text,k,profile):
    maxProb = -1.;
    for kmer in getSubKmers(text,k):
        p = getProfProb(kmer,profile,k)
        if p > maxProb:
            maxProb = p
            pmp_kmer = kmer
    return pmp_kmer

def printGPMPK(filename):
    with open(filename)as f:
        text = f.readline().strip()
        k = int(f.readline().strip())
        profile =[ [0. for j_col in range(k)] for i_row in range(4)]
        for i_row in range(4):
            profile[i_row] = [float(el) for el in
                              f.readline().strip().split(' ')]
    print(getProfileMostProbKmer(text,k,profile))

def getMotifScore(motifs):
    score = 0
    t = len(motifs)
    k = len(motifs[0])
    for j_col in range(k):
        count = [0 for i in range(4)]
        for i_row in range(t):
            count[ref.find(motifs[i_row][j_col])] += 1
        score += t - max(count)
    return score  

def getProfileFromMotif(motifs,pseudoCount=0):
    t = len(motifs)
    k = len(motifs[0])
    profile = [[pseudoCount for i in range(k)] for j in range(4)]
    for i_row in range(t):
        for j_col in range(k):
            profile[ref.find(motifs[i_row][j_col])][j_col]+=1
    for j_col in range(k):
        tot = float(sum(profile[i][j_col] for i in range(4)))
        for i in range(4):
            profile[i][j_col] /= tot
    return profile

def getGreedySearchMotif(dna,k,t):
    bestMotifs = [next(getSubKmers(dna[i],k)) for i in range(t)]
    bestScore = getMotifScore(bestMotifs)
    for kmer in getSubKmers(dna[0],k):
        motifs = [kmer]
        for text in [dna[i] for i in range(1,t)]:
            profile = getProfileFromMotif(motifs,1)
            motifs.append(getProfileMostProbKmer(text,k,profile))
        score = getMotifScore(motifs)
        if score < bestScore:
            bestScore = score
            bestMotifs = motifs
    return bestMotifs

def printGGSM(filename):
    with open(filename) as f:
        [k,t] = [int(el) for el in f.readline().strip().split(' ')]
        dna = [el.strip() for el in f.readlines()]
    for el in getGreedySearchMotif(dna,k,t):
        print(el)

from random import randint
def getRamdomSearchMotif(dna,k,t,n_iter=1000):
    bestOuterScore = 1e10
    for i_iter in xrange(n_iter):
        ik =  [randint(0,len(dna[0])-k) for i in range(t)]
        motifs = [dna[i][ik[i]:ik[i]+k] for i in range(t)]
        bestInnerMotifs = motifs
        bestInnerScore = getMotifScore(motifs)
        while True:
            profile = getProfileFromMotif(motifs,1)
            motifs = [getProfileMostProbKmer(text,k,profile) for text in dna]
            score = getMotifScore(motifs)
            if score >= bestInnerScore:
                break
            bestInnerScore = score
            bestInnerMotifs = motifs
        if bestInnerScore < bestOuterScore:
            bestOuterScore = bestInnerScore
            bestMotifs = bestInnerMotifs          
    return bestMotifs
        
def printGRSM(filename):
    with open(filename) as f:
        [k,t] = [int(el) for el in f.readline().strip().split(' ')]
        dna = [el.strip() for el in f.readlines()]
    for el in getRamdomSearchMotif(dna,k,t):
        print(el )   

from random import random, randint
def getProbWeightedRandKmer(text,k,profile):
    kmers = [el for el in getSubKmers(text,k)]
    p = [getProfProb(kmer,profile,k) for kmer in kmers]
    ptot = sum(p)
    cumuP = [(p[i] + sum(p[0:i]))/ptot for i in range(len(p))]
    pick=random()
    return kmers[[el>pick for el in cumuP].index(True)]


def getGibbsSampMotif(dna,k,t,Ninner=100,Nouter=20):
    bestOuterScore = 1e10
    for i_out in xrange(Nouter):
        ik =  [randint(0,len(dna[0])-k) for i in range(t)]
        motifs = [dna[i][ik[i]:ik[i]+k] for i in range(t)]
        bestInnerMotifs = motifs[:]
        bestInnerScore = getMotifScore(motifs)
        for i_in in xrange(Ninner):
            i = randint(0,t-1)
            profile = getProfileFromMotif([motifs[j] for j in xrange(t)
                                           if j != i],1)
            motifs[i] = getProbWeightedRandKmer(dna[i],k,profile)
            score = getMotifScore(motifs)
            if score < bestInnerScore:                
                bestInnerScore = score
                bestInnerMotifs = motifs[:]
        if bestInnerScore < bestOuterScore:
            bestOuterScore = bestInnerScore
            bestMotifs = bestInnerMotifs[:]
    return bestMotifs

def printGGSM(filename):
    with open(filename) as f:
        [k,t,n_in] = [int(el) for el in f.readline().strip().split(' ')]
        dna = [el.strip() for el in f.readlines()]
    M = getGibbsSampMotif(dna,k,t,n_in)
    for el in M:
        print(el)
    return M


##        for j ← 1 to N
##            i ← Random(t)
##            Profile ← profile matrix constructed from all strings in Motifs
##                       except for Motifi
##            Motifi ← Profile-randomly generated k-mer in the i-th sequence
##            if Score(Motifs) < Score(BestMotifs)
##                BestMotifs ← Motifs
##        return BestMotifs

#GGGCCGTTGGT GGACCGTTGAC  
##cInp = readClumpInp('dataset_4_5.txt')
##print ' '.join([i for i in
##                getClumpKmers(cInp['Genome'],cInp['k'],cInp['t'],cInp['L'])])

##with open('E-coli.txt') as f:
##    clumps = getClumpKmers(f.readline().strip(),9,3,500)

##sampleFile = 'dataset_2_6.txt'
##with open(sampleFile) as f:
##    [imp,pat] = [l.strip() for l in f.readlines()]
##print("The pattern {0:s} appears in the target {1} times.".format(
##        pat,PatternCount(imp,pat)))
##sampleFile = 'dataset_2_9.txt'
##with open(sampleFile) as f:
##    [imp,k] = [l.strip() for l in f.readlines()]
##k = int(k);
##
##print(FrequentWords(imp,int(k)))
##print(FreqWords(imp,int(k)))

##sampleFile = 'dataset_3_2.txt'
##with open(sampleFile) as f:
##    imp = f.readline().strip()
##print(revComp(imp))

##sampleFile = 'dataset_3_5.txt'
##sampleFile = 'Vibrio_cholerae.txt'
##with open(sampleFile) as f:
##    imp = f.readline().strip()
##    print(findOcc(imp, 'CTTGATCAT'))

    

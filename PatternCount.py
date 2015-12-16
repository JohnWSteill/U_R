def PatternCount(Text, Pattern):
    count = srchWindowBegin = 0        
    while True:
        srchWindowBegin = Text.find(Pattern, srchWindowBegin) + 1
        if srchWindowBegin > 0:
            count+=1
        else:
            return count

def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = []
    for i in range(0,len(Text)-k):
        Pattern = Text[i:i+k];
        Count.append( PatternCount(Text, Pattern))
    maxCount = max(Count)
    for i in range(0,len(Text)-k):
        if (Count[i]==maxCount and not (Text[i:i+k] in FrequentPatterns)):
           FrequentPatterns.append(Text[i:i+k])
    return FrequentPatterns

ref = 'ACGT'
def strToKmerInd(kmer):
    return sum([ref.find(kmer[i])*(4**i) for i in range(k)])

def kmerIndToStr(val):
    retStr = ['-' for i in range(k)]
    for i in range(k):
        [val, retStr[i]] = [val/4, ref[val%4]]
    return ''.join(retStr)

def FreqWords(Text,k): # Using Index
    count = [ 0 for i in range(4**k)]
    for i in range(len(Text)-k):
        count[strToKmerInd(Text[i:i+k])]+=1
    maxCount = max(count)
    inds = [ ind for ind in range(len(count)) if (count[ind] == maxCount)]
    FrequentPatterns = [kmerIndToStr(i) for i in inds]
    return FrequentPatterns
    
def revComp(Text):
    comp = dict([('A','T'), ('T','A'), ('C','G'), ('G','C')])
    return ''.join([comp[z] for z in Text[::-1]])
    
def findOcc(Text, Pattern):
    inds = []
    srchWindowBegin = 0        
    while True:
        srchWindowBegin = Text.find(Pattern, srchWindowBegin) + 1
        if srchWindowBegin > 0:
            inds.append(srchWindowBegin-1)
        else:
            return ' '.join(str(i) for i in inds)

def getClumpKmers(Genome,k,t,L):
    clumpKmers = []
    clumpInds = [0 for i in range(4**k)]
    
    

##BetterClumpFinding(Genome, k, t, L)
##        FrequentPatterns ← an empty set
##        for i ←0 to 4k − 1
##            Clump(i) ← 0
##        Text ← Genome(0, L)
##        FrequencyArray ← ComputingFrequencies(Text, k)
##        for i ← 0 to 4k − 1
##            if FrequencyArray(i) ≥ t
##                Clump(i) ← 1
##        for i ← 1 to |Genome| − L
##            FirstPattern ← Genome(i − 1, k)
##            index ← PatternToNumber(FirstPattern)
##            FrequencyArray(index) ← FrequencyArray(index) − 1
##            LastPattern ← Genome(i + L − k, k)
##            index ← PatternToNumber(LastPattern)
##            FrequencyArray(index) ← FrequencyArray(index) + 1
##            if FrequencyArray(index) ≥ t
##                Clump(index) ← 1
##        for i ← 0 to 4k − 1
##            if Clump(i) = 1
##                Pattern ← NumberToPattern(i, k)
##                add Pattern to the set FrequentPatterns
##        return FrequentPatterns


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

sampleFile = 'dataset_3_5.txt'
sampleFile = 'Vibrio_cholerae.txt'
with open(sampleFile) as f:
    imp = f.readline().strip()
    print(findOcc(imp, 'CTTGATCAT'))

    

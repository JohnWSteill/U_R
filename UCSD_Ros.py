# http://docs.sublimetext.info/en/latest/basic_concepts.html
class UCSD_Ros_Solver(object): 
    """ This class performs Rosalind Project solutions. 

    Usage is as follows:  
    r = UCSD_Ros_Solver() 
    r.myProb("MyData.txt")
    Assumes that MyData.txt is in working directory or getwd()/Data 
    """      

    def StringRecon(self,inpFile,echoToScreen=False):
        """ Solve the String Composition Problem.

        Input: An integer k and a string Text.
        Output: Compositionk(Text) (the k-mers can be provided in any order).

        Sample Input:
        5
        CAATCCAAC

        Sample Output:
        CAATC
        AATCC
        ATCCA
        TCCAA
        CCAAC
        """
        outFileName = inpFile + "_out"
        with open(inpFile) as f:
            k = int(f.readline().strip())
            dna = f.readline().strip()
        
        kMers = set()

        for substr in [dna[i:i+k] for i in xrange(len(dna)-k+1)]:
            kMers.add(substr)

        with open(outFileName,"w") as f:
            for kMer in kMers:
                f.write(kMer + "\n")     



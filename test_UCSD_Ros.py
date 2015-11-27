#UCSD_Ros_test.py
import unittest
import os
import tempfile
import filecmp

import UCSD_Ros
reload(UCSD_Ros)

def makeTmpFile(self,data,fName):
	with open(fName,"w") as f:
		for d in data:
			f.write(d)
	with open(fName) as f:
		for line in f:
			print ("ins: ", line)

class testLoader(unittest.TestCase):

	def test1(self):
		testInput = "5 \nCAATCCAAC" 
		testOutput = "CAATC\nAATCC\nATCCA\nTCCAA\nCCAAC"
		makeTmpFile(self,testInput,"tmpTestInp.txt")
		makeTmpFile(self,testOutput,"tmpTestOut.txt")
		rs = UCSD_Ros.UCSD_Ros_Solver()
		rs.UCSD_StringRecon("tmpTestInp.txt")
		self.assertTrue(filecmp.cmp("tmpTestOut.txt","tmpTestInp.txt"+"_out"))

	def test2(self):
		self.assertTrue(True)


if __name__ == "__main__":
	unittest.main()


##print 'Building a file name yourself:'
##filename = '/tmp/guess_my_name.%s.txt' % os.getpid()
##temp = open(filename, 'w+b')
##try:
##    print 'temp:', temp
##    print 'temp.name:', temp.name
##finally:
##    temp.close()
##    # Clean up the temporary file yourself
##    os.remove(filename)
##
##print
##print 'TemporaryFile:'
##temp = tempfile.TemporaryFile()
##try:
##    print 'temp:', temp
##    print 'temp.name:', temp.name
##finally:
##    # Automatically cleans up the file
##    temp.close()
#sclass testLoader
# import os
# import tempfile

# directory_name = tempfile.mkdtemp()
# print directory_name
# # Clean up the directory yourself
# os.removedirs(directory_name)

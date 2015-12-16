#UCSD_Ros_test.py
import unittest
import os
import tempfile
import filecmp
import glob
import UCSD_Ros
reload(UCSD_Ros)

def makeTmpFile(data,fName):
	with open(fName,"w") as f:
		for d in data:
			f.write(d)

class testLoader(unittest.TestCase):

	def test1(self):
		for (testInput, testOutput) in [
				("5 \nCAATCCAAC" ,"CAATC\nAATCC\nATCCA\nTCCAA\nCCAAC\n"),
				("5 \nCAATC" ,"CAATC\n")]:			
			makeTmpFile(testInput,"tmpTestInp.txt")
			makeTmpFile(testOutput,"tmpTestOut.txt")
			rs = UCSD_Ros.UCSD_Ros_Solver()
			rs.StringRecon("tmpTestInp.txt")
			with open("tmpTestOut.txt") as f:
				setTruth = set([el.strip() for el in f.readlines()])
			with open("tmpTestInp.txt"+"_out") as f:
				setTest = set([el.strip() for el in f.readlines()])
			self.assertTrue(setTest==setTruth)
			print(setTruth,setTest)
			for tempf in glob.glob('./tmp*'):
				os.remove(tempf)


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

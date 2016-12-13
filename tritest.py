from triarray import *
import numpy as np


testArr = TriArr(6)

for x in [1,18,3,4,5,6]:
	testArr.addValue(x)
print testArr.toString()
for x in [6,5,4,6,2,5,17,3,1,3,5,4,1,-8,3]:
	testArr.addValue(x)
print testArr.segmentTri(4,0)
print testArr.toString()
for z in xrange(12):
	x,y=testArr.linTo2D(z)
	print z,x,y,testArr.getValue(x,y),testArr.data_array[z]
#exit()
print testArr.data_array
print testArr.toString()

subTriangle = testArr.getSubTriangle(1,4)
for z in xrange(6):
	x, y = subTriangle.linTo2D(z)
	print z,x,y,subTriangle.getValue(x,y),subTriangle.data_array[z]

print subTriangle.data_array
print subTriangle.toString()

subTriangle.data_array = [1,2,3]
print subTriangle.toString()



for pos in xrange(0,12):
	print pos,testArr.data_array[pos], testArr.linTo2D(pos)

print '\n\n'
region=np.array([1,2,3,4])
tri_arr = TriArr(region.shape[0])
for x in xrange(region.shape[0]):
	for y in xrange(x, region.shape[0]):
		val = x+y*0.1#np.sum(region[x:y]) / np.sqrt(y - x)
		tri_arr.addValue(val)
		print x,y,val,tri_arr.getValue(x,y)
print tri_arr.data_array
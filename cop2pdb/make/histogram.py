import sys







#argv1 is input files argv2 is number of bins argv3 is output number
f = open(sys.argv[1], 'r')
fl = f.readlines()
fls = [line.split() for line in fl]
flsx = [float(line[0]) for line in fls]
#strips down the file, and makes a list of all the values
f.close()

binnum = int(sys.argv[2])
valmin = 0
valmax = float(sys.argv[4])

binsize = (valmax-valmin) / float(binnum)


hist = []
for i in range(binnum):
  hist.append(0)
  for j in range(len(flsx)):
    if flsx[j] >= binsize*i and flsx[j] < binsize*(i+1):
      hist[i] += 1
for i in flsx:
  if float(i) == valmax:
    print('yo')
    hist[-1] += 1
if sum(hist) != 252:
  print("sum not the same as data points")
  print(sys.argv[1])
  print(sum(hist))
g = open("hist" + sys.argv[3], 'w')
for i in range(len(hist)):
  g.write(str(i*binsize) + " " + str(hist[i]) + "\n")
g.close()
  



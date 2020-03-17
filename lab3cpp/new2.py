from matplotlib import pyplot as plt
x = []
y = []
xx = []
yy = []
with open('Newton_two.dat', encoding='utf-8') as file:
 for line in file:
        x1, y1 = line.split()
        x.append([float(x1)])
        y.append([float(y1)])
file.close
with open('orignew2.dat', encoding='utf-8') as file:
 for line in file:
        x2, y2 = line.split()
        xx.append([float(x2)])
        yy.append([float(y2)])
file.close
fig1 = plt.figure()
fig2 = plt.figure()
axt = fig1.add_subplot()
axt.plot(yy,xx,color="blue", marker = "*")
axt.plot(y,x,color="green")
fig1.set_figwidth(6)
fig1.set_figheight(6)

axm = fig2.add_subplot()
axm.plot(yy,xx,color="red" , marker="*")
fig1.set_figwidth(6)
fig1.set_figheight(6)

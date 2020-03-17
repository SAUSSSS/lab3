import numpy as np
import matplotlib.pyplot as plt
import csv
import math

def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i+=step
        
def a_cof(u,n):
    temp = u
    for i in range(n):
        temp *= (u - i)
    return temp

def a_co(u,n):
    temp = u
    for i in range(n):
        temp *= (u + i)
    return temp

def fact(n):
    f = 1;
    for i in range(2, n + 1):
        f *= i;
    return f; 



row = []
tablet = []
def csv_reader(file_obj):

    reader = csv.reader(file_obj)
    for row in reader:
            row.append(" ".join(row))
            tablet.append(row)
if __name__ == "__main__":
    csv_path = "dat.csv"
    with open(csv_path, "r") as f_obj:
        csv_reader(f_obj) 
  


def lagrang(x,y,dx):
    result = 0.0
    for i in range(len(x)):
        P = 1.0
        for j in range(len(x)):
            if(j != i):
                P *= (dx - y[j]) / (y[i] - y[j])
        if(x[i] != 999.9):
            result += P * x[i]
    if(result < 100):
        return result
   
    
def newton_two(x,y,dx): 
    n = len(x)
    mas = [[0 for i in range(len(x))] for j in range(len(x))]
    for i in range(len(x)):
        mas[i][0] = x[i]
    j = n
    for i in range(1,len(x)):
        for t in range(len(x)):
             if j > i:
                 mas[t][i] = mas[t][i - 1] - mas[t - 1][i - 1];
        j -= 1
    sum = mas[n - 1][0]
    u = (dx - y[n - 1]) / (y[1] - y[0])
    for i in range(1,len(x)):
        sum += (a_co(u,i) * mas[n - 1][i]) / fact(i)
    if(sum < 100):
        return sum 
    
    
 
def newton_one(x,y,dx): 
    mas = [[0 for i in range(len(x))] for j in range(len(x))]
    for i in range(len(x)):
        mas[i][0] = x[i]
    for i in range(1,len(x)):
        for j in range(len(x) - i):
             mas[j][i] = mas[j + 1][i - 1] - mas[j][i - 1]; 
    sum = mas[0][0]
    u = (dx - y[0]) / (y[1] - y[0])
    for i in range(1,len(x)):
        sum += (a_cof(u,i) * mas[0][i]) / fact(i)
    if(sum < 100):
        return sum
    



 
x = []
y = []
z = []
res = []
print("Enter amount of data to analyze: ")
size = float(input())
print("row and column") 
rrow = int(input())
col = int(input())



for i in range(int(size)):
    x.append(float(tablet[i+rrow][col]))
    y.append(i)
    z.append(i)
    
#тут надо менять вызов функции(lagrange,newton_one,newton_two)    
res.append([newton_two(x,y,i) for i in frange(0,size,0.1)])
yy = []
t = 0.0
for i in range(len(res[0])):
    yy.append(t)
    t += 0.1
resn = []   

#resn.append([*reversed(res[0])]) #Для newton_one  
    

print(resn) 
fig1 = plt.figure()
fig2 = plt.figure()
axt = fig1.add_subplot()
axz = fig2.add_subplot()
axt.plot(yy,res[0],color="blue")
axt.plot(y,x,color="green")
fig1.set_figwidth(6)
fig1.set_figheight(6)



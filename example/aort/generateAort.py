import numpy
f=open('aort.bdy','w')
start=float(input("Введите начало аорты отсчета на сетке ->\n"))
end  =float(input("Ввежите конец  аорты на сетке->\n"))
h    =float(input("Введите высоту аорты->\n"))
y1   =float(input("Введите координату нижней границы аорты->\n"))
y2   =float(input("Введите координату верхней границы аорты->\n"))
N    =float(input("Введите количество точек->\n"))
array = numpy.linspace(int(start),int(end),int(N))
writeCoords = lambda a,up: [f.write(str(x) + '\t' + str(up) + '\n') for x in a]
writeCoords(array,y1)
writeCoords(array,y1 + h)
f.close()

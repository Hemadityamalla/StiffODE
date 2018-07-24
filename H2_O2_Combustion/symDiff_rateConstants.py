from sympy import symbols, diff, init_printing

init_printing()

ko,ki = symbols("k,K")

y1, y2, y3, y4, y5, y6, y7, y8 = symbols("a,b,c,d,e,f,g,h")

Pr = (ko/ki)*(2.5*y1 + y2 + 16*y3 + y4 + y5 + y6 + y7 + y8)

f = (ko*y2*y4)/(1 + Pr)
print('Symbol notation:\n')
print("ko=k \n ki=K \n y1=a \n y2=b \n y3=c \n y4=d \n y5=e \n y6=f \n y7=g \n y8=h \n \n")



print('-----------Given function---------------- \n',f, '\n')

print("-------------------Differentiated:----------\n")

diff(f, y1)

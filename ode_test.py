import scipy as s
import scipy.integrate as intg
t_inc=0.001


def integral(t):
	""""""
	global t_inc
	return(intg.odeint(funcion,1,[t,t+t_inc]))
	
def funcion(t,y):
	""""""
	return y	 	

def main():
	""""""
	global t_inc
	resultados=s.zeros([2,6000])
	i=0
	for t in arange(-1,5,t_inc):
		resultados[0][i] = t
		[kk, resultados[1][i]] = integral(t)
		i=i+1
	return resultados
	
if __name__ == '__main__':
	res = main()

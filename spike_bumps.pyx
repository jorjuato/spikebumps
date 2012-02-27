import scipy as s
import pickle as pickle


							#####################################################
t = 0.000					# Tiempo de simulacion								#
t_sim = 1.0					# Tiempo total de duración de la simulación			#
t_inc = 0.001				# Incremento temporal								#
							#####################################################
							
							#####################################################
N = 100						# Tamaño de la red									#
P = 3						# Número de patrones almacenados					#
							#####################################################


#################################################################################
#			GLOSARIO DE VARIABLES EMPLEADAS EN EL PROGRAMA						#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
################################################################################



def iniciar_experimento():
	""""""
	import scipy as s
	exp = Experimento()
	exp.__train__()
	#exp.__test__()
	file = open('experimento','w')
	pi = pickle.Pickler(file)
	pi.dump(exp)
	file.close()
	return exp



class Experimento:
	"""
		Clase que dirige el desarrollo de las fases de entrenamiento y testeo
		de la red. Contiene los métodos y datos nesarios para la extraccion y 
		análisis de los resultados

	"""
	
	def __init__(self):	
	
		#global t 								
		#global t_sim 					
		#global t_inc 					
		#global N						
		#global P	
		self.minCue=0.2
		self.maxCue=1
		self.intCue=0.2		
										#################################################
		self.N = N 						# Número total de neuronas en la red
		self.P = P   					# Número de patrones almacenados
		self.a = 0.2					# Memory sparseness
		self.M = 10						# Factor normalizador de la plasticidad sináptica
										##################################################
										
										################################################
		self.cueMax = 1					# Fidelidad del patron presentado
		self.cueMin = 0.2				#
		self.cueInc = 0.4				#
		self.kMax = 23					# Numero medio de conexiones salientes por neurona
		self.kMin = 3					#
		self.kInc = 10					#
		self.qMax = 1					# Parametro de aleatoriedad de las conexiones
		self.qMin = 0					#
		self.qInc = 0.2					###############################################
		
		self.lambda_0 = 0.25
		self.lambda_cue = 0.1


										###############################################
		self.t0 = 0.15					# Cue onset time
		self.t1 = 0.3					# Tiempo de inicio del decrecimiento del cue
		self.t2 = 0.5					# Cue offset time
		self.t_win = 50 				# Sampling time window		
		self.t_sim = t_sim				################################################
		self.t_inc = t_inc
		
		####### Generacion aleatoria de patrones ###############################
		self.patrones = s.random.rand(self.P, self.N) 
		for i in range(self.P) :					  
			temp = s.absolute(s.random.rand(self.N) - (0.5 - self.a)) #Usar binomial
			self.patrones[i] = temp.round()
		########################################################################
		self.red = Red(self.N)
		self.resultados = []
		self.k=0
		self.q=0
		self.parametros = Parametros()
		self.parametros.cueMax = self.cueMax
		self.parametros.cueMin = self.cueMin
		self.parametros.cueInc = self.cueInc
		self.parametros.kMax = self.kMax
		self.parametros.kMin = self.kMin
		self.parametros.kInc = self.kInc
		self.parametros.qMax = self.qMax
		self.parametros.qMin = self.qMin
		self.parametros.qInc = self.qInc
		self.parametros.tMax = t_sim
		self.parametros.tInc = t_inc
		self.parametros.N = N
		self.parametros.P = P
		self.parametros.patrones = self.patrones[:]
		
		#### Declaracion de la lista de resultados  #################################
		for k in range(int((self.kMax-self.kMin+self.kInc)/self.kInc)):
			self.resultados.append([])
			for q in range(int((self.qMax-self.qMin+self.qInc)/self.qInc)):
				self.resultados[k].append(Trial())
				for mu in range(self.P):
					for cue in range(int((self.cueMax-self.cueMin+self.cueInc)/self.cueInc)):
						#for neu in self.red.neuronas :
						self.resultados[k][q].patrones[mu].test_rates.append([])
						self.resultados[k][q].patrones[mu].test_spikes.append([])
		##############################################################################
						
	def secuencia_simulaciones(self):
		""

		for k in range(self.kMin,self.kMax,self.kInc):
			for q in s.arange(self.qMin,self.qMax,self.qInc):
				# Generar la red, conectarla e inicializar variables 
				self.red = Red(self.N)
				self.__connect__(self.q)
				self.red.edta = s.zeros((self.N,self.P))
				self.red.Edta_media = 0.000
				self.red.J = s.zeros((self.N,self.N))	
				# Entrenamos la red generada
				self.__train__()
				self.__edta1__()
				self.__J1__()				
				# Testamos para cada patron
				self.__test__()
				# Analizamos los datos obtenidos
				self.__analize__()
				#Estos dos incrementos finales permiten llevar un contador
				#global a la clase del punto en que se encuentra la simulacion.
				self.q = self.q + 1
			self.k = self.k + 1
				
	def __connect__(self, q):
		"""
			Funcion que determina las conexiones de la red a través de una ley
			de probalidad gobernada por el parámetro q. Segun su valor, la red 
			oscila entre aleatoria y simétrica
		"""
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				if ni.id == nf.id: continue
				self.red.conexiones[ni.id,nf.id]= round(s.random.random()) # binomial y tal
	
		
	def __train__(self):
		""
		for mu in  range(self.P):
			self.__train_results__(mu, self.__show_pattern__(self.patrones[mu]))

			
			
	def __test__(self):
		"Fase de testeo del algoritmo."	
		for mu in range(self.P) :
			for cue in arange(self.cueMin,self.cueMax,self.cueInc):				
				self.__test_results__(mu, cue, self.__show_pattern__(self.__build_pattern__(self.patrones[mu]), cue))
	
	def __analize__(self):
		""
		
	def __build_pattern__(self, _patron):
		""	
		inicio = int(s.random.random_integers(self.N))
		fin = inicio + int(cueQ*self.N)
		if fin >= self.N :
			for i in range(self.N-fin,inicio) :
				_patron[i] = 0
		else :
			for i in range(0,inicio) :
				_patron[i] = 0
			for i in range(fin, self.N) :
				_patron[i] = 0
		return _patron	

	def __show_pattern__(self, _patron):
		"""
			Para todos los patrones se repite la siguiente secuencia de entrenamiento:
				1) Aleatorizamos los V y dejamos simulando 200 ms
				2) Imponemos el patrón durante 200ms
				3) Decaemos el patron durante 200ms
				4) Simulamos hasta completar el tiempo total 
		"""
		spikes_temp = []
		for neu in self.red.neuronas :
			neu.V = neu.V_res + (neu.V_thr-neu.V_res)/(s.random.rand() * 10)
		for self.t in s.arange(0, self.t0, t_inc) :
			spikes_temp.append(self.__simula__(s.zeros(self.N) ))
		for self.t in s.arange(self.t0, self.t1, t_inc) :
			spikes_temp.append(self.__simula__(_patron))
		for self.t in s.arange(self.t1, self.t2, t_inc) :
			_patron = _patron* \
			(self.lambda_cue * (self.t2-self.t1)/(self.t2-self.t1)) + self.lambda_0
			spikes_temp.append(self.__simula__(_patron))
		_patron = s.zeros(self.N)	
		for self.t in s.arange(self.t2, self.t_sim, t_inc) :
			spikes_temp.append(self.__simula__(_patron))
		return spikes_temp
			
	def __simula__(self, _patron):	
		""
		registro = []
		ns = self.red.neuronas
		for neu in ns :
			registro.append(neu.simula(_patron[neu.id], self.red.I_syn, self.red.r_t))
		return registro
		
	def __edta1__(self, mu):	
		""
		
		
	def __edta2__(self):	
		""
		self.edta

	def __J1__(self, mu):
		"Calcula la matriz de interacciones hebbianas"
		res = self.resultados[self.k][self.q]
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				for index in range(0, mu) :
					res.J[ni.id,nf.id] = res.J[ni.id,nf.id] + \
					(res.edta[ni.id,index]/res.edta_media - 1) * \
					(res.edta[nf.id,index]/res.edta_media - 1)
		self.J = self.J / M	

	def __J2__(self):
		"Calcula la matriz de interacciones hebbianas"
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				for index in range(0, self.P) :
					self.J[ni.id,nf.id] = self.J[ni.id,nf.id] + \
					(self.edta[ni.id,index]/self.edta_media - 1) * \
					(self.edta[nf.id,index]/self.edta_media - 1)
		self.J = self.J / M	
		
	
			
	def __train_results__(self, mu, spikes):
		""
		result = self.resultados[self.k][self.q]
		result.patrones[mu].conexiones = self.red.conexiones[:]
		result.J = self.red.J[:]
		result.patrones[mu].edta = self.red.edta[:]
		result.patrones[mu].edta_media = self.red.edta_media[:]
		
		rates = result.patrones[mu].train_rates
		spike_times = result.patrones[mu].train_spikes
		inc = int(self.t_sim/self.t_win)
		
		for neu in self.red.neuronas : 
			rates.append([])
			spike_times.append(neu.spikecounts[:])
			for cont in range(inc):
				rates.append(sum(float(spikes[neu.id][cont:cont+inc]))/inc)


	def __test_results__(self, mu, cue, spikes):
		""
		
		spike_times = self.resultados[self.k][self.q].patrones[mu].test_spikes
		inc = self.t_sim/self.t_win
		for cue in range(int((self.cueMax-self.cueMin+self.cueInc)/self.cueInc)):
			rates = self.resultados[self.k][self.q].patrones[mu].test_rates[cue]
			for neu in self.red.neuronas :
				rates.append([])
				spike_times[cue].append(neu.spikecounts[:])
				for cont in range(inc):
					rates[neu.id].append(sum(spikes[neu.id][cont:cont+inc]))
				
		
class Red:
	"""
		Clase que encapsula el conjunto de neuronas y variables globales
		necesarios para definir el comportamiento de la red	
	"""

	def __init__(self, N):
		#global t					# Tiempo de simulacion					
		#global t_sim 				# Tiempo total de duración de la simulación
		#global t_inc 				# Incremento temporal
		
		self.neuronas = []
		for i in range(N):
			self.neuronas.append(Neurona(i))
		self.conexiones = s.zeros((N,N))
		#self.input = s.zeros(N)
		self.I_syn = s.zeros(N)
		self.r_t = 0
		self.N = N
		self.lambda_syn = 40
		self.J = s.zeros((N,N))
		self.edta = s.zeros((P,N))
		self.edta_media = s.zeros(P)
		self.tau_1 = 0.03
		self.tau_2 = 0.004
	
	def __I_syn__(self):
		"Calcula una matriz de interacciones sinapticas"
		self.I_syn=s.zeros(self.N)
		for nf in self.neuronas:
			for ni in self.neuronas:
				if self.conexiones[ni.id,nf.id] == 0 : continue
				for spike in ni.spikecounts:
					if t - spike > 0.2 : break				
					self.I_syn[nf.id] = self.I_syn[nf.id] + \
								tau_m/(R_m*(tau_1-tau_2))*\
								(s.exp(-(t-spike)/tau_1) - s.exp(-(t-spike)/tau_2))*\
								self.lambda_syn * (ni.V_thr - ni.V_res)/N *\
								(1+self.J(ni.id,nf.id)) 
			
	def __r_t__(self): # cuidadin con cuando se llama a esta funcion
		"Calcula la actividad media de la red"
		r_i, r_media = 0, 0
		for neu in self.neuronas:
			r_media = r_media + r_i
			r_i = 0
			for spike in neu.spikecounts[::-1]:				
				if t - spike > t_inc*2 : 
					break	
				else :	
					r_i += r_i
		return float(r_media) / self.N
		



class Neurona :
	""""Unidad básica del modelo"""
	
	def __init__(self,i):
		#global t 					# Tiempo de simulacion					
		#global t_sim					# Tiempo total de duración de la simulación
		#global t_inc				# Incremento temporal
		
		self.id = i  # Contiene el índice que nos permitirá reconocer a la neurona.
		self.V = 0.00
		self.I_tot = 0.00
		self.I_inh = 0.00
		self.I_t = 0.00
		self.spikecounts = []
		
		#self.lambda_syn = 40
		self.lambda_inh = 20
		self.lambda_0 = 0.25
		self.lambda_cue = 0.1
		
		self.tau_m = 0.005

		self.tau_inh = 0.004
		self.tau_ref = 0.003
		
		self.V_thr = 1.00
		self.V_res = 0.00
		self.R_mem = 1.00
		
	def simula(self, _patron, I_syn, r_t):
		"Lleva a cabo un paso de la simulacion"
		self.I_tot = self.__I_inh__(r_t) + self.__I_t__(_patron) + I_syn 
		self.V = self.__V_t__()
		if len(self.spikecounts) != 0 :
			if self.__refractario__() : 
				return 0
			elif self.V >= self.V_thr :
				self.spikecounts.append(t)
				self.V = self.V_res
				return 1
			else : return 0

	def __refractario__(self):
		"Comprueba que la ultima spiga no se prudujo hace menos de t_ref"
		if (t - self.spikecounts[-1]) < self.t_ref : 	
			return 1
		else: 
			return 0	
		
	cdef __I_inh__(self, r_t):
		
		return s.exp(-t/self.tau_m)*r_t	#Hay que calcular bien esa integral	

	cdef __I_t__(self, _patron):
		"Lleva a cabo un paso de la simulacion"
		return (self.lambda_0 + self.lambda_cue) * (self.V_thr - self.V_res)


	cdef __V_t__(self):
		
		return self.V + s.exp(self.I_tot*t)




class Parametros:
	pass
		

class Trial: 
	def __init__(self):
		global P
		self.J = []
		self.conexiones = [[]]
		self.edta = []
		self.edta_media = []
		self.patrones = []
		for mu in range(P):
			self.patrones.append(Patron())
			self.edta.append([])
		
class Patron:
	def __init__(self):
		self.train_rates = []
		self.train_spikes = []
		self.test_rates = []
		self.test_spikes = []

if __name__ == '__main__':
	iniciar_experimento()

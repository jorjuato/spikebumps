

import scipy as sci

#cdef public void initspikes()

t = 0.000					# Tiempo de simulacion					
t_sim = 1					# Tiempo total de duración de la simulación
inc_t = 0.001				# Incremento temporal

N = 50 					# Número total de neuronas en la red
k = 24   					# Numero medio de conexiones salientes por neurona
q = 0    					# Parametro de aleatoriedad de las conexiones

P = 3						# Número de patrones almacenados
a = 0.2						# Memory sparseness
M = 10						# Factor normalizador de la plasticidad sináptica

t0 = 0.15					# Cue onset time
t1 = 0.3					# Tiempo de inicio del decrecimiento del cue
t2 = 0.5					# Cue offset time
t_win = 50 					# Sampling time window

lambda_syn = 40
lambda_inh = 20
lambda_0 = 0.25
lambda_cue = 0.1
minCue=0.2
maxCue=1
intCue=0.2
tau_m = 0.005
tau_1 = 0.03
tau_2 = 0.004
tau_inh = 0.004
tau_ref = 0.003

V_thr = 1.00
V_res = 0.00
R_mem = 1.00

def iniciar_experimento():
	""
	exp = Experimento()
	exp.train()
	exp.test()
	return exp
	

	
# ###########################################################################
cdef class Neurona:
	"Unidad básica del modelo"
	
	def __init__(self,i):
		self.id = i  # Contiene el índice que nos permitirá reconocer a la neurona.
		self.V = 0.00
		self.I_tot = 0.00
		self.I_inh = 0.00
		self.I_t = 0.00
		self.spikecounts = []

	def simula(self, _patron, I_syn, r_t):
		"Lleva a cabo un paso de la simulacion"
		self.I_inh = self.__I_inh__(r_t)
		self.I_t = self.__I_t__(_patron)
		self.I_tot = self.I_inh + self.I_t + I_syn 
		self.V = self.__V_t__()
		if len(self.spikecounts) != 0 :
			if self.__refractario__() : 
				return
			elif self.V >= V_thr :
				self.spikecounts.append(t)
				self.V = V_res

	def __refractario__(self):
		"Comprueba que la ultima spiga no se prudujo hace menos de t_ref"
		if (t - self.spikecounts[-1]) < t_ref : 	
			return 1
		else: 
			return 0	
		
	def __I_inh__(self, r_t):
		"Lleva a cabo un paso de la simulacion"
		return exp(-t/tau_m)*r_t	#Hay que calcular bien esa integral	

	def __I_t__(self, _patron):
		"Lleva a cabo un paso de la simulacion"
		return (lambda_0 + lambda_cue) * (V_thr - V_res)


	def __V_t__(self):
		"Retorna el valor del voltaje en el instante de tiempo"
		return self.V + exp(self.I_tot*t)


# ###########################################################################

cdef class Red:
	"""
		Clase que encapsula el conjunto de neuronas y variables globales
		necesarios para definir el comportamiento de la red	
	"""

	def __init__(self, N):

		self.neuronas = []
		for i in range(N):
			self.neuronas.append(Neurona(i))
		self.conexiones = zeros((N,N))
		
		self.Input = zeros(N)
		self.I_syn = zeros(N)
		self.r_t = 0
		self.registros = [] 
		
		
	def __I_syn__(self):
		"Calcula una matriz de interacciones sinapticas"
		self.I_syn=zeros(N)
		for nf in self.neuronas:
			for ni in self.neuronas:
				if self.conexiones[ni.id,nf.id] == 0 : continue
				for spike in ni.spikecounts:
					if t - spike > 1 : break				
					self.I_syn[nf.id] = I_syn[nf.id] + \
								tau_m/(R_m*(tau_1-tau_2))*\
								(exp(-(t-spike)/tau_1) - exp(-(t-spike)/tau_2))*\
								lambda_syn * (V_thr - V_res)/N *\
								(1+self.J(ni.id,nf.id)) 
			
	def __r_t__(self): # cuidadin con cuando se llama a esta funcion
		"Calcula la actividad media de la red"
		r_i, r_media = 0, 0
		for neu in self.neuronas:
			r_media = r_media + r_i
			r_i = 0
			for spike in neu.spikecounts[::-1]:				
				if t - spike > inc_t*2 : 
					break	
				else :	
						r_i = r_i + 1
		return r_media
		
	def simula(self, patron):
		"""
		"""
		for neu in self.neuronas :
			neu.simula(patron[neu.id], self.I_syn, self.r_t)
			#self.Edta[neu.id, patron[neu.id]] = 0####Como calcularlo??
			
		#self.Edta_media = mean(self.Edta)
		#__J__()
	

	def record(self):
		"Maneja las estructuras de datos necesarias para guardar los resultados"
		for neu in self.neuronas:
			rate = 0			
			for spike in neu.spikecounts[::-1] :
				if t - spike > t_win: 
					break
				else:
					rate = rate + 1
			self.registros[neu.id].append(rate)
			#([t neu.V ,neu.I_tot ,neu.I_syn ,neu.I_inh ,neu.I_t ,neu.r_t])



# ###########################################################################
cdef class Experimento:
	"""
		Clase que dirige el desarrollo de las fases de entrenamiento y testeo
		de la red. Contiene los métodos y datos nesarios para la extraccion y 
		análisis de los resultados

	"""
	
	def __init__(self):	
		self.patrones = sci.random.rand(P, N)
		for i in range(P) :
			temp = abs(sci.random.rand(N) - (0.5 - a)) #Usar binomial
			self.patrones[i] = temp.round()
		self.red = Red(N)
		self.__connect__(q)
		self.edta = zeros((N,P))
		self.Edta_media = 0.000
		self.resultados = [[]]
		self.J = zeros((N,N))		
		
	def train(self):
		"""
			Fase de entrenamiento del algoritmo. Para todos los patrones se repi
			te la siguiente secuencia de entrenamiento:
				1) Aleatorizamos los V y dejamos simulando 200 ms
				2) Imponemos el patrón durante 200ms
				3) Decaemos el patron durante 200ms
				4) Simulamos hasta completar el tiempo total 
		"""
		for mu in  range(P):
			self.__show_pattern__(self.patrones[mu])
			self.__edta__(mu)
			self._J(mu)
			
	def test(self):
		"""
			Fase de testeo del algoritmo. 
		"""
		for cueQ in sci.arange(minCue,maxCue,intCue) :
			for mu in range(P) :
				_patron = self.__build_pattern__(cueQ, self.patrones[mu])
				self.__show_pattern__(_patron)
				self.resultados.append(self.__obtain_results__(mu, cueQ))
				
				
	def __connect__(self, q):
		"""
			Funcion que determina las conexiones de la red a través de una ley
			de probalidad gobernada por el parámetro q. Segun su valor, la red 
			oscila entre aleatoria y simétrica
		"""
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				if ni.id == nf.id: continue
				self.red.conexiones[ni.id,nf.id]= round(sci.random.random()) # binomial y tal
	
	def __edta__(self, mu):	
		""
		self.edta

	def _J(self, mu):
		"Calcula la matriz de interacciones hebbianas"
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				for index in range(0, mu) :
					self.J[ni.id,nf.id] = self.J[ni.id,nf.id] #+ \
					#(self.edta[ni.id,index]/self.edta_media - 1) * \
					#(self.edta[nf.id,index]/self.edta_media - 1)
		self.J = self.J / M		
		
	def __show_pattern__(self, _patron):
		""
		
		#print "Entrenando el patron numero : ",_patron 
		for neu in self.red.neuronas :
			neu.V = V_res + (V_thr-V_res)/(sci.random.rand() * 10)
		for t in sci.arange(0, 0.2, inc_t) :
			self.red.simula(zeros(N))
		#	print "iteracion: ", t , "\n"
		for t in sci.arange(t, t+0.2, inc_t) :
			self.red.simula(_patron)
		#	print "iteracion: ", t , "\n"
		for t in sci.arange(t, t+0.2, inc_t) :
			_patron = _patron*(lambda_cue*(t-t1)/(t2-t1)) + lambda_0
			self.red.simula(_patron)
		#	print "iteracion: ", t , "\n"
		_patron = zeros(N)	
		for t in sci.arange(t, t_sim, inc_t) :
			self.red.simula(_patron)
		#	print "iteracion: ", t , "\n"			
	
	def __obtain_results__(self, mu, cueQ):
		""
		pass
		
	def __build_pattern__(self, cueQ, _patron):
		""	
		rnd = int(sci.random.random_integers(N))
		tot = rnd + int(cueQ*N)
		if tot >= N :
			for i in range(tot-N,rnd) :
				_patron[i] = 0
		else :
			for i in range(0,rnd) :
				_patron[i] = 0
			for i in range(tot, N) :
				_patron[i] = 0
		return _patron
		
		
# ###########################################################################

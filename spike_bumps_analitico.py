import scipy as s
import pickle as pickle


#################################################################################
#						VARIABLES GLOBALES										#
#################################################################################
t = 0.000					# Tiempo 											#
t_sim = 1.0					# Tiempo total de duración de la simulación			#
t_inc = 0.001				# Incremento temporal								#
N = 2						# Tamaño de la red									#
P = 1						# Número de patrones almacenados					#
#################################################################################


#################################################################################
#					PARAMETROS FISIOLOGICOS DE LA SIMULACION					#
#################################################################################
lambda_syn = 40				# constante 										#
lambda_inh = 20				# constante 										#
lambda_0 = 0.025			# constante 										#
lambda_cue = 0.1			# constante 										#
#################################################################################
tau_inh = 0.004				# constante temporal de inhibicion					#
tau_ref = 0.003				# constante temporal del periodo refractario		#
#################################################################################													#		
V_thr = -50.00				# Potencial umbral									#
V_res = -70.00				# Potencial de membrana en reposo					#
R_mem = 1.00				# Resistencia de la membrana						#
tau_m = 0.005				# constante de tiempo de la membrana				#
#################################################################################
tau_1 = 0.03				# constante temporal sinaptica rapida				#	
tau_2 = 0.004				# constante temporal sinaptica lenta				#	
#################################################################################



def iniciar_experimento():
	""""""
	import scipy as s
	exp = Experimento()
	exp.secuencia_simulaciones()
	#del exp.red	
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
	
	###################  VARIABLES GLOBALES #############################################
		self.N = N 						# Número total de neuronas en la red			#
		self.P = P   					# Número de patrones almacenados				#
		self.t_sim = t_sim				# Tiempo total de simulacion					#
		self.t_inc = t_inc				# Incremento temporal							#
	#####################################################################################
							
							
	###############  PARAMETROS GENERALES DEL EXPERIMENTO  ##############################
		self.cueMax = 1					# Fidelidad del patron presentado				#	
		self.cueMin = 0.2				#												#
		self.cueInc = 0.4				#												#
		self.kMax = 23					# Numero medio de conexiones por neurona		#
		self.kMin = 3					#												#
		self.kInc = 20					#												#
		self.qMax = 1					# Parametro de aleatoriedad de las conexiones	#
		self.qMin = 0					#												#
		self.qInc = 0.5					#												#
		self.a = 0.2					# Memory sparseness								#
		self.M = 10						# Factor normalizador de plasticidad sináptica	#
	#####################################################################################
		
		
	######################### 	CUE PARAMETERS  #########################################
		self.t0 = 0.15					# Cue onset time								#
		self.t1 = 0.3					# Tiempo de inicio del decrecimiento del cue	#
		self.t2 = 0.5					# Cue offset time								#
		self.t_win = 0.05				# Sampling time window							#
		self.lambda_0 = 0.25			# Corriente basal de entrada					#
		self.lambda_cue = 0.1			# Factor aplicado a la corriente del cue		#
	#####################################################################################
	
	
	########### GENERACION ALEATORIA DE PATRONES ########################################
		self.patrones = s.random.rand(self.P, self.N) 									#
		for i in range(self.P) :					  									#
			temp = s.absolute(s.random.rand(self.N) - (0.5 - self.a)) #Usar binomial	#
			self.patrones[i] = temp.round()												#
	#####################################################################################
	
	
	#########  INSTANCIACION DE CLASES FUNDAMENTALES DEL EXPERIMENTO  ###################	
		self.red = Red(self.N)															#
		self.resultados = []															#
		self.parametros = Parametros()													#
		self.k=0																		#
		self.q=0																		#
	#####################################################################################
	
	
	############ DECLARACION DEL LISTADO DE RESULTADOS  #################################
		for k in range(int((self.kMax-self.kMin+2*self.kInc)/self.kInc)):				#
			self.resultados.append([])													#
			for q in range(int((self.qMax-self.qMin+2*self.qInc)/self.qInc)):			#
				self.resultados[k].append(Trial())										#
				for mu in range(self.P):												#
					for cue in range(int((self.cueMax-self.cueMin+self.cueInc)/self.cueInc)):		#
						self.resultados[k][q].patrones[mu].test_rates.append([])		#
						self.resultados[k][q].patrones[mu].test_spikes.append([])		#
	#####################################################################################
	

	#################	ALMACENADO LOCAL DE PARAMETROS  #################################		
		self.parametros.cueMax = self.cueMax											#
		self.parametros.cueMin = self.cueMin											#
		self.parametros.cueInc = self.cueInc											#
		self.parametros.kMax = self.kMax												#
		self.parametros.kMin = self.kMin												#
		self.parametros.kInc = self.kInc												#
		self.parametros.qMax = self.qMax												#
		self.parametros.qMin = self.qMin												#
		self.parametros.qInc = self.qInc												#
		self.parametros.tMax = t_sim													#
		self.parametros.tInc = t_inc													#
		self.parametros.N = N															#
		self.parametros.P = P															#
		self.parametros.patrones = self.patrones[:]										#
	#####################################################################################

						
	def secuencia_simulaciones(self):
		"""
			Funcion que se encarga de llamar secuencialmente a todas las funciones 
			necesarias para llevar a cabo el experimento. Lleva el control de los
			parámetros k y q.
		"""
		for k in range(self.kMin,self.kMax+self.kInc,self.kInc):
			self.q = 0
			for q in s.arange(self.qMin,self.qMax+self.qInc,self.qInc):
				# Conectamos las neuronas de la red
				self.__connect__(k,q)
				
				# Entrenamos la red generada
				self.__train__()
				self.__edta1__()
				self.__J1__()				
				
				# Testamos para cada patron
				self.__test__()
				
				# Analizamos los datos obtenidos
				self.__analize__()
				
				# Estos dos incrementos finales permiten llevar un contador
				# global a la clase del punto en que se encuentra la simulacion.
				self.q = self.q + 1
			self.k = self.k + 1
				
	def __connect__(self, k=None, q=None):
		"""
			Funcion que determina las conexiones de la red a través de una ley
			de probalidad gobernada por el parámetro q. Segun su valor, la red 
			oscila entre aleatoria y simétrica
		"""
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				if ni.id == nf.id: continue
				self.red.conexiones[ni.id,nf.id] = int(round(s.random.random() - self.a)) 
				# binomial y tal
				# usar k y q
		
	def __train__(self):
		"Fase de entrenamiento del algoritmo"
		for mu in  range(self.P):
			self.__train_results__(mu, self.__show_pattern__(self.patrones[mu]))

	def __test__(self):
		"Fase de testeo del algoritmo."	
		for mu in range(self.P) :
			for cue in s.arange(self.cueMin,self.cueMax+self.cueInc,self.cueInc):				
				self.__test_results__(mu, cue, self.__show_pattern__(self.__build_pattern__(self.patrones[mu], cue)))
	
	def __analize__(self):
		"Analisis de la respuesta de la red ante las distintas condiciones experimentales"
		
	def __build_pattern__(self, _patron, cueQ):
		"Funcion encargada de generar un patron con una fidelidad determinada por cueQ"	
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
		global t
		spikes_temp = []
		for neu in self.red.neuronas :
			neu.V = neu.V_res + (neu.V_thr-neu.V_res)/(s.random.rand() * 10)
		for t in s.arange(0, self.t0, t_inc) :
			spikes_temp.append(self.__simula__(s.zeros(self.N) ))
		for t in s.arange(self.t0, self.t1, t_inc) :
			spikes_temp.append(self.__simula__(_patron))
		for t in s.arange(self.t1, self.t2, t_inc) :
			#_patron = _patron* \
			#(self.lambda_cue * (self.t2-self.t1)/(self.t2-self.t1)) + self.lambda_0
			spikes_temp.append(self.__simula__(_patron))
		_patron = s.zeros(self.N)	
		for t in s.arange(self.t2, self.t_sim, t_inc) :
			spikes_temp.append(self.__simula__(_patron))
		return spikes_temp
			
	def __simula__(self, _patron):	
		"""
			Para un instante de tiempo, esta funcion actualiza los valores del sistema
			simulando el paso de inc_t segundos. Devuelve la actividad de la red como
			una lista de unos y ceros
		"""
		registro = s.zeros(N,s.float_)
		ns = self.red.neuronas
		self.red.__I_syn__()
		self.red.r_t = self.red.__r_t__()
		for neu in ns :
			registro[neu.id]=neu.simula(_patron[neu.id], self.red.I_syn[neu.id], self.red.r_t)
		return registro
		
	def __edta1__(self):	
		"""
			Calcula las frecuencias de disparo tanto individuales como medias para cada
			patrón pero cuando todos han sido mostrados a la red
		"""
		res = self.resultados[self.k][self.q]
		for mu in range(self.P):
			pat = res.patrones[mu]
			for neu in self.red.neuronas:
				# Cogemos ultimos (intervalo/t_win) rates y hacemos la media
				last = -int(0.150/self.t_win)
				lista = pat.train_rates[neu.id][last :]
				res.edta[mu,neu.id] = sum(lista[:])
			res.edta_media[mu] = sum(res.edta[mu,:])
			
		
	def __edta2__(self, mu0):	
		"""
			Calcula las frecuencias de disparo tanto individuales como medias para cada
			patrón actualizandose cada vez que se enseña un nuevo patron.
		"""
		res = self.resultados[self.k][self.q]
		for mu in range(0, mu0):
			pat = res.patrones[mu]
			for neu in self.red.neuronas:
				# Cogemos ultimos (intervalo/t_win) rates y hacemos la media
				res.edta[mu,neu.id] = sum(pat.train_rates[-int(0.150/self.t_win) :])
			res.edta_media[mu] = sum(res.edta[mu,:])
			
	def __J1__(self):
		"""
			Calcula la matriz de interacciones hebbianas despues de enseñar todos los 
			patrones a la red
		"""
		res = self.resultados[self.k][self.q]
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				for mu in range(self.P) :
					res.J[ni.id,nf.id] = res.J[ni.id,nf.id] + \
					(res.edta[mu, ni.id]/res.edta_media[mu] - 1) * \
					(res.edta[mu, nf.id]/res.edta_media[mu] - 1)
		res.J = res.J / self.M	

	def __J2__(self, mu0):
		"""
			Calcula la matriz de interacciones hebbianas despues de enseñar todos los 
			patrones a la red
		"""
		res = self.resultados[self.k][self.q]
		for ni in self.red.neuronas :
			for nf in self.red.neuronas :
				for mu in range(0, mu0) :
					res.J[ni.id,nf.id] = res.J[ni.id,nf.id] + \
					(res.edta[mu, ni.id]/res.edta_media[mu] - 1) * \
					(res.edta[mu, nf.id]/res.edta_media[mu] - 1)
		res.J = res.J / self.M	
		
	
			
	def __train_results__(self, mu, spikes):
		"""
			Calcula y almacena los resultados de la fase de entrenamiento. Toma como 
			parametros el numero de patron y la actividad booleana de la red
		"""
		result = self.resultados[self.k][self.q]
		result.conexiones = self.red.conexiones.copy()
	
		rates = result.patrones[mu].train_rates
		spike_times = result.patrones[mu].train_spikes
		inc = int(self.t_sim*self.t_win)
		
		for neu in self.red.neuronas : 
			rates.append([])
			spike_times.append(neu.spikecounts[:])
			for cont in range(inc):
				rates.append(sum(float(spikes[neu.id][cont:cont+inc]))/inc)


	def __test_results__(self, mu, cue, spikes):
		"""
			Calcula y almacena los resultados de la fase de testeo para cada patron. 
			Toma como parametros el numero de patron, la calidad del estimulo
			y la actividad booleana de la red
		"""
		inc = int(self.t_sim*self.t_win) 			# Numero de ventanas en la simulacion
		
		for cue in range(int((self.cueMax-self.cueMin+self.cueInc)/self.cueInc)):
			rates = self.resultados[self.k][self.q].patrones[mu].test_rates[cue]
			spike_times = self.resultados[self.k][self.q].patrones[mu].test_spikes[cue]
			for neu in self.red.neuronas :
				rates.append([])
				spike_times.append(neu.spikecounts[:])
				for cont in range(inc):
					rates[neu.id].append(sum(spikes[neu.id][cont:cont+inc]))
				
		
class Red:
	"""
		Clase que encapsula el conjunto de neuronas y variables globales
		necesarios para definir el comportamiento de la red	
	"""
	global t
	
	def __init__(self, N):
	
	################# GENARACION DE LA LISTA DE NEURONAS Y SUS CONEXIONES  ##############
		self.neuronas = []																#
		for i in range(N):																#
			self.neuronas.append(Neurona(i))											#
		self.conexiones = s.zeros((N,N))												#
		self.r_t = 0																	#
		self.I_syn = s.zeros(N, s.float_)
		self.J = s.zeros([N,N], s.float_)		#
	#####################################################################################	
	
	################# Almacenamiento local de paramatros fisiologicos	#################
		self.tau_1 = tau_1																#
		self.tau_2 = tau_2																#
		self.lambda_syn = lambda_syn													#
		self.V_thr = V_thr																#
		self.V_res = V_res 																#
		self.R_mem = R_mem																#
		self.tau_m = tau_m																#
	#####################################################################################
	
	def __I_syn__(self):
		"Calcula un vector de interacciones sinapticas"
		global t
		self.I_syn=self.I_syn*0
		const = self.lambda_syn*(self.V_thr - self.V_res)/(N*(self.tau_1-self.tau_2))
		for nf in self.neuronas:
			for ni in self.neuronas:
				if self.conexiones[ni.id,nf.id] == 0 : continue
				for spike in ni.spikecounts:
					if t - spike > 0.2 : break				
					self.I_syn[nf.id] = self.I_syn[nf.id] +\
								(s.exp(-(t-spike)/self.tau_1) - s.exp(-(t-spike)/self.tau_2))*\
								(1+self.J[ni.id,nf.id])
		self.I_syn = self.I_syn*(self.tau_m/self.R_mem)
		
	def __r_t__(self): 
		"Calcula la actividad media de la red"
		global t
		r_i, r_media = 0, 0
		for neu in self.neuronas:
			r_media = r_media + r_i
			r_i = 0
			for spike in neu.spikecounts[::-1]:				
				if t - spike > t_inc*2 : 
					break	
				else :	
					r_i += r_i
		return (float(r_media)) / N
		



class Neurona :
	""""Unidad básica del modelo"""
	
	global t
	
	def __init__(self,i):
		self.id = i  # Contiene el índice que nos permitirá reconocer a la neurona.
		self.V = V_res
		self.voltaje = []
		self.I_tot = 0.00
		#self.I_inh = 0.00
		#self.I_t = 0.00
		self.spikecounts = []
		
		self.lambda_inh = lambda_inh
		self.lambda_0 = lambda_0
		self.lambda_cue = lambda_cue
		
		self.tau_inh = tau_inh
		self.tau_ref = tau_ref
		
		self.V_thr = V_thr
		self.V_res = V_res
		self.R_mem = R_mem
		self.tau_m = tau_m
		
	def simula(self, _patron, I_syn, r_t):
		"Lleva a cabo un paso de la simulacion"
		global t
		self.I_tot = -self.__I_inh__(r_t) + self.__I_t__(_patron) + I_syn 
		self.voltaje.append(self.V)
		self.V = self.__V_t__()
		if self.__refractario__() == 1: 
			return 0
		if self.V >= self.V_thr :
			self.spikecounts.append(t)
			self.V = self.V_res
			return 1
		return -1
	
	def __refractario__(self):
		"Comprueba que la ultima spiga no se prudujo hace menos de t_ref"
		global t
		if len(self.spikecounts) > 0 :
			if (t - self.spikecounts[-1]) < self.tau_ref : 	
				self.V = self.V_res
				return 1
			else: 
				return 0
		else:
			return 4				
		
	def __I_inh__(self, r_t):
		"Lleva a cabo un paso de la simulacion"
		global t
		return s.exp(-t/self.tau_m)*r_t	#Hay que calcular bien esa integral	

	def __I_t__(self, _patron):
		"Lleva a cabo un paso de la simulacion"
		global t
		return (self.lambda_0 + self.lambda_cue) * (self.V_thr - self.V_res)


	def __V_t__(self):
		"Retorna el valor del voltaje en el instante de tiempo"
		global t
		return self.V + s.exp(self.I_tot*t*self.R_mem/self.tau_m)*0.0001



class Parametros:
	pass
		

class Trial: 
	def __init__(self):
		self.J = s.zeros((N,N), s.float_)
		self.conexiones = s.zeros((N,N))
		self.edta = s.zeros((P,N), s.float_)
		self.edta_media = s.zeros((P), s.float_)
		self.patrones = []
		for mu in range(P):
			self.patrones.append(Patron())
				
class Patron:
	def __init__(self):
		self.train_rates = []
		self.train_spikes = []
		self.test_rates = []
		self.test_spikes = []

if __name__ == '__main__':
	experimento = iniciar_experimento()

l = 44
i = 12345 
class MiClase:
	"Simple clase de ejemplo"
	global i
	def __init__(self):
		self.k = "kakaka"	
		
	def f(self):
		
		i = i + 1
		print i , l, self.k

import time
import functools

def elapsedTime(string : str):
	def chrono(func) :
		@functools.wraps(func)
		def wrapper(*args, **kwargs) :
			start = time.process_time_ns()
			valeur = func(*args, **kwargs)
			stop = time.process_time_ns()
			elapsed = (stop - start) / 10**9
			print(string)
			print(elapsed)
			return valeur
		return wrapper
	return chrono
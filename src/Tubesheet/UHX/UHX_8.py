
from ._UHX_common import Configuration

class UHX_8:
	def __init__(self, shellsideflange, tubesideflange):
		self._W_c = tubesideflange.W
		self._W_m_1_c = tubesideflange.W_m_1
		self._W_m_1_s = shellsideflange.W
		self._W_s = shellsideflange.W_m_1
		self._W_m_1_max = max(self._W_m_1_c, self._W_m_1_s)
		self._W_max = max(self._W_c, self._W_s)
		self.load_case_choices = ["Design 1", "Design 2", "Design 3", "Design 4","Operating 1", "Operating 2", "Operating 3", "Operating 4"]

		self._table = {
			Configuration.a: {"Design 1": 0, 			 "Design 2": 0, 			"Design 3": 0, 				 "Design 4": 0, "Operating": 0}			 ,
			Configuration.b: {"Design 1": self._W_m_1_c, "Design 2": 0, 			"Design 3": self._W_m_1_c,   "Design 4": 0, "Operating": self._W_c}  ,
			Configuration.c: {"Design 1": self._W_m_1_c, "Design 2": 0, 			"Design 3": self._W_m_1_c,   "Design 4": 0, "Operating": self._W_c}  ,
			Configuration.d: {"Design 1": self._W_m_1_c, "Design 2": self._W_m_1_s, "Design 3": self._W_m_1_max, "Design 4": 0, "Operating": self._W_max},
			Configuration.e: {"Design 1": 0, 			 "Design 2": self._W_m_1_s, "Design 3": self._W_m_1_s, 	 "Design 4": 0, "Operating": self._W_s}	 ,
			Configuration.f: {"Design 1": 0, 			 "Design 2": self._W_m_1_s, "Design 3": self._W_m_1_s,   "Design 4": 0, "Operating": self._W_s}	 ,
			Configuration.A: {"Design 1": 0, 			 "Design 2": 0, 			"Design 3": 0, 				 "Design 4": 0, "Operating": 0}			 ,
			Configuration.B: {"Design 1": self._W_m_1_c, "Design 2": 0, 			"Design 3": self._W_m_1_c, 	 "Design 4": 0, "Operating": self._W_c}  ,
			Configuration.C: {"Design 1": self._W_m_1_c, "Design 2": 0, 			"Design 3": self._W_m_1_c, 	 "Design 4": 0, "Operating": self._W_c}  ,
			Configuration.D: {"Design 1": 0, 			 "Design 2": 0, 			"Design 3": 0, 				 "Design 4": 0, "Operating": 0}
		}

	def _w_star(self, configuration, case):
		try:
			return self._table[configuration][case]
		except KeyError:
			return None

import numpy as np
import matplotlib.pyplot as plt



class ParksMcClellan():
	def __init__(self, L, wp, ws, delta1, delta2, epsilon=1e-3):
		self.L = L
		self.wp = wp
		self.ws = ws
		self.del1 = delta1
		self.del2 = delta2
		self.epsilon = epsilon
		self.npoints = int(1e3)
		self.w = np.arange(0, np.pi, 1e-3)
		self.alter_w = None
		self.Wwk = None
		self.opt_del = None
		self.Hd = None
		self.Ae = None
		self.nitr = 0

	def hd(self, omega):
		Hd = np.ones(omega.shape)
		Hd[omega>self.wp] = 0
		return Hd

	def weight_wk(self):
		self.Wwk = np.zeros(self.alter_w.shape)
		self.Wwk[self.alter_w >= self.ws] = 1
		self.Wwk[self.alter_w <= self.wp] = self.del2/self.del1

	def b(self):
		B = np.zeros(shape=self.alter_w.shape)
		for i in range(len(B)):
			bi = 1
			for j in range(len(B)):
				if i!=j:
					bi *= 1/(np.cos(self.alter_w[j])-np.cos(self.alter_w[i]))
			B[i] = bi
		return B

	def Opt_del(self):
		bk = self.b()
		a = np.dot(bk, self.Hd)
		b  = 0
		for i in range(len(bk)):
			b += bk[i]*np.power(-1, i+1)/self.Wwk[i]
		return a/b

	def C(self):
		self.Hd = self.hd(self.alter_w)
		self.weight_wk()
		alter_ones = np.power(-1, np.arange(self.Hd.shape[0]))
		self.opt_del = self.Opt_del()
		C = self.Hd - alter_ones*self.opt_del/self.Wwk
		return C

	def d(self):
		return (self.b()*(np.cos(self.alter_w)-np.cos(self.alter_w[-1])))[:-1]

	def Ae_ejw(self):
		a = 0
		b = 0
		d = self.d()
		C = self.C()
		A = np.zeros(self.w.shape)
		#for j in range(w.shape[0]):
		for i in range(self.alter_w.shape[0]-1):
			# if(np.absolute(w[j]-self.alter_w[i]) <= 1e-3):
			# 	a = C[i]
			# 	b = 1
			# 	break
			# else:
			a+=d[i]*C[i]/(np.cos(self.w)-np.cos(self.alter_w[i]))
			b+=d[i]/(np.cos(self.w)-np.cos(self.alter_w[i]))
		A = a/b
		# alter_ones = np.power(-1, np.arange(self.Hd.shape[0]))
		# print(alter_ones)
		# extreme_case = alter_ones.copy().astype('float')
		# for i in range(self.L+2):
		# 	if i < int(self.L/2)+2:
		# 		extreme_case[i] = 1 - (self.del1/self.del2)*self.opt_del*alter_ones[i]
		# 	else:
		# 		extreme_case[i] = self.opt_del*alter_ones[i]
		# print(extreme_case)
		# print((self.del1/self.del2)*self.opt_del)
		# idxs = []
		# rounded_w = np.round(self.w, 3)
		# for i in np.round(self.alter_w, 3):
		# 	idxs.append(np.where(rounded_w==i)[0][0])
		# A[idxs] = extreme_case
		# for i in range(1, A.shape[0]-1):
		# 	if (A[i+1] - A[i])/A[i] >= 2:
		# 		A[i+1] = A[i]
		return A

	def gd(self, y, x, init, type=0):
		# Type 1--> local Maxima, 0--> local Minima
		x_ = init
		init_idx = np.where(x==init)[0]

	def find_new_set(self):
		W = np.zeros(self.w.shape)
		for j in range(W.shape[0]):
			if self.w[j] <= self.wp:
				W[j] = self.del2/self.del1
			elif self.w[j] >= self.ws:
				W[j] = 1
		E = W*(self.hd(self.w)-self.Ae)
		#req_space = np.logical_or(self.w<=self.wp, self.w>=self.ws)
		#E = E[req_space]
		return E#, self.w[req_space]


	def iter(self, init_guess=None, show=False):

		if (init_guess==None).all():
			Nin_wp = int(self.L/2)+1
			Nin_ws = self.L + 2 -Nin_wp
			#np.random.seed(1500)
			self.alter_w = np.append(np.random.choice(self.w[self.w<=self.wp][1:], Nin_wp), np.random.choice(self.w[self.w>=self.ws][:-1], Nin_ws))
			self.alter_w = np.sort(self.alter_w)
		else:
			self.alter_w = init_guess
		self.Ae = self.Ae_ejw()
		print(self.C())
		alter_ones = np.power(-1, np.arange(self.Hd.shape[0]))
		# M = np.power(np.cos(self.alter_w[0]), np.arange(self.L+1))
		# for i in range(1, self.L+2):
		# 	M = np.vstack((M, np.power(np.cos(self.alter_w[i]), np.arange(self.L+1))))
		# M = np.hstack((M, (alter_ones/self.Wwk).reshape(self.L+2, 1)))
		# x = np.linalg.inv(M)@(self.Hd.reshape(self.L+2, 1))
		# print(self.opt_del, x[-1])
		print(self.Wwk)
		self.nitr+=1
		e = self.find_new_set()
		# print((self.del2/self.del1)*self.opt_del)
		# extreme_case = alter_ones
		# for i in range(self.L+2):
		# 	if i < int(self.L/2)+2:
		# 		extreme_case[i] = 1 + (self.del2/self.del1)*self.opt_del*alter_ones[i]
		# 	else:
		# 		extreme_case[i] = self.opt_del*alter_ones[i]
		# print(extreme_case)
		# idxs = []
		# rounded_w = np.round(self.w, 3)
		# for i in np.round(self.alter_w, 3):
		# 	idxs.append(np.where(rounded_w==i)[0][0])
		if show:
			plt.figure()
			plt.title(str(self.nitr)+"th Iteration")
			plt.plot(self.w, self.Ae)
			aws = self.w[self.w>=self.ws]
			bwp = self.w[self.w<=self.wp]
			plt.plot(bwp, 1-self.opt_del*np.ones(bwp.shape), "r--")
			plt.plot(bwp, 1+self.opt_del*np.ones(bwp.shape), "r--")
			plt.plot(bwp, 1+self.del1*np.ones(bwp.shape), "b--")
			plt.plot(bwp, 1-self.del1*np.ones(bwp.shape), "b--")
			plt.plot(aws, -1*self.opt_del*np.ones(aws.shape), "r--")
			plt.plot(aws, self.opt_del*np.ones(aws.shape), "r--")
			plt.plot(aws, self.del2*np.ones(aws.shape), "b--")
			plt.plot(aws, -1*self.del2*np.ones(aws.shape), "b--")
			plt.figure()
			plt.plot(bwp, self.del1*np.ones(bwp.shape), "r--")
			plt.plot(bwp, -1*self.del1*np.ones(bwp.shape), "r--")
			plt.plot(aws, self.del2*np.ones(aws.shape), "r--")
			plt.plot(aws, -1*self.del2*np.ones(aws.shape), "r--")
			plt.plot(self.w, e)
			#plt.plot(self.alter_w, self.Ae[idxs], "g*")
			plt.show()

pi = np.pi
init_guess = np.append(np.linspace(1e-3, pi/2, 6), np.linspace(2*pi/3, pi, 6, endpoint=False))
#init_guess = np.array([0, pi/12, pi/6, pi/4, pi/3, 5*pi/12, 2*pi/3, 2*pi/3+pi/18, 2*pi/3+pi/9, 2*pi/3+pi/6, 2*pi/3+2*pi/9, 2*pi/3+5*pi/18])
h_n = ParksMcClellan(10, np.pi/2, 2*np.pi/3, delta1=0.0125, delta2=0.025)
h_n.iter(init_guess=init_guess, show=True)
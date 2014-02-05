 #!/usr/bin/python
 # -*- coding: utf-8 -*-
 
##### Authorship ##############################
###############################################
# @name	The X-mu Python Library
# @author	Daniel Lewis
# @email	daniel@vanirsystems.com
# @website	http://vanirsystems.com/
# @version	alpha-4
# @license	Apache License, Version 2.0
""" Copyright 2014 Daniel Lewis
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
###############################################

##### Requirements ############################
from sympy import *
import sympy.mpmath as mpmath
from sympy.assumptions.assume import *
from numpy.random import rand
###############################################
	
ALPHA = Symbol('alpha', positive=True, real=True, bounded=True)
global_assumptions.add(Q.is_true(And((ALPHA >= 0.0), (ALPHA <= 1.0))))

class Xmu(object):
	""" The Xmu class (which extends the Python object class), is a somewhat abstract class that sets up an object for polymorphism, and also provides the basic properties (such as x, u, mu and the Xmu function) """
	
	x = Symbol('x', real=True, bounded=True)
	muequals = Piecewise((0.0, True))
	xequals = EmptySet()
	u = Interval(0.0, 1.0)
	
	def __init__(self, u):
		""" __init__ contructor taking in a singule parameter u representing the universe. Sets the object u, and also sends assertions to SymPy.
		@param	u	an interval being a subset of R (real numbers), representing the universe of this particular domain.
		@return	An instantiated X-mu object """
		
		self.set_u(u)
		Q.is_true((self.x >= self.u.start) & (self.x <= self.u.end))
	
	def __str__(self):
		""" __str__ method returns a string representation of the X-mu formula. 
		@return	string """
		
		return str(self.get_muequals())
	
	def setMuFunction(self, func):
		""" Sets the mu function. Currently a synonym for set_muequals.
		@param	func	a sympy function """
		self.set_muequals(func)
	
	def set_muequals(self, func):
		""" Sets the mu function.
		@param	func	a sympy function """
		self.muequals = func
	
	def get_muequals(self):
		""" Gets the mu function.
		@return	muequals	a sympy function """
		return self.muequals
		
	def setXFunction(self, func):
		""" Sets the X-mu function. Currently a synonym for set_xequals.
		@param	func	a sympy function """
		self.set_xequals(func)
		
	def set_xequals(self, func):
		""" Sets the X-mu function.
		@param	func	a sympy function """
		self.xequals = func
	
	def get_xequals(self):
		""" Gets the X-mu function.
		@return	xequals	a sympy function """
		return self.xequals
	
	def set_u(self, u):
		""" Sets the universe u.
		@param	u	a sympy interval """
		self.u = u
	
	def intersectX(self, target):
		""" Performs X-mu set intersection (itself union a target), and returns the result.
		@param	target	an Xmu instance. 
		@return	intersection	an BasicXmu instance """
		i1 = target.get_xequals()
		i2 = self.get_xequals()
		return BasicXmu(self.u, i1.intersect(i2))
	
	def differenceX(self, target):
		""" Performs X-mu set difference (itself - a target), and returns the result.
		@param	target	an Xmu instance. 
		@return	difference	an BasicXmu instance """
		i1 = target.get_xequals()
		i2 = self.get_xequals()
		result = i2 - i1
		return BasicXmu(self.u, result)
	
	def unionX(self, target):
		""" Performs X-mu set union (itself union a target), and returns the result.
		@param	target	an Xmu instance. 
		@return	union	an BasicXmu instance """
		i1 = target.get_xequals()
		i2 = self.get_xequals()
		return BasicXmu(self.u, i1.union(i2))
	
	def negateX(self):
		""" Performs X-mu set negation (- itself), and returns the result.
		@todo	Currently not working. This is due to SymPy not quite working as expected in terms of bounded intervals and symbolics therein.
		@return	negation	an BasicXmu instance """
		result = self.u - self.get_xequals()
		return BasicXmu(self.u, result)
	
	def arithmeticalOperationX(self, i2, alpha, operation):
		""" Wrapper class for performing Fuzzy Arithmetic on X-mu Functions using SymPy/MPI. Primarily used as a private class, but could be used publicly.
		@param	i2	the target sympy formula.
		@param	alpha	the level on which to perform the operation.
		@param	operation	One of the following: +, -, *, /, **
		"""
		res = self.get_xequals().subs(ALPHA, float(alpha))
		res2 = i2.get_xequals().subs(ALPHA, float(alpha))
		if (not res.is_EmptySet) and (not res2.is_EmptySet):
			if (type(res) != FiniteSet) and (type(res2) != FiniteSet):
				result = eval("res.to_mpi() " + operation + " res2.to_mpi()")
			else:
				res_mpi = mpmath.mpi(res.inf, res.sup)
				res2_mpi = mpmath.mpi(res2.inf, res2.sup)
				result = eval("res_mpi " + operation + " res2_mpi")
			
			return BasicXmu(self.u, Interval(float(mpmath.mpf(result.a)), float(mpmath.mpf(result.b))))
		else:
			return None
	
	def multiplyX(self, i2, alpha):
		""" Performs X-mu Multiplication at a certain alpha point.
		@param	i2	a target interval 
		@param	alpha	the specified alpha point
		@note	multiplyX() is a non-symbolic function due to SymPy interactions. Therefore, it must be done at a specific alpha point. """
		return self.arithmeticalOperationX(i2, alpha, "*")
	
	def powX(self, i2, alpha):
		""" Performs X-mu Power at a certain alpha point.
		@param	i2	a target interval 
		@param	alpha	the specified alpha point
		@note	powX() is a non-symbolic function due to SymPy interactions. Therefore, it must be done at a specific alpha point. """
		return self.arithmeticalOperationX(i2, alpha, "**")
	
	def addX(self, i2, alpha):
		""" Performs X-mu Addition at a certain alpha point.
		@param	i2	a target interval 
		@param	alpha	the specified alpha point
		@note	addX() is a non-symbolic function due to SymPy interactions. Therefore, it must be done at a specific alpha point. """
		return self.arithmeticalOperationX(i2, alpha, "+")
	
	def subX(self, i2, alpha):
		""" Performs X-mu Subtraction at a certain alpha point.
		@param	i2	a target interval 
		@param	alpha	the specified alpha point
		@note	subX() is a non-symbolic function due to SymPy interactions. Therefore, it must be done at a specific alpha point. """
		return self.arithmeticalOperationX(i2, alpha, "-")
	
	def divX(self, i2, alpha):
		""" Performs X-mu Division at a certain alpha point.
		@param	i2	a target interval 
		@param	alpha	the specified alpha point
		@note	divX() is a non-symbolic function due to SymPy interactions. Therefore, it must be done at a specific alpha point. """
		return self.arithmeticalOperationX(i2, alpha, "/")
		
class BasicXmu(Xmu):
	""" BasicXmu is a wrapper class around an Xmu function so that we can perform set theoretic operations, thereupon. All X-mu Set Theoretic operations return an instance of this BasicXmu class. """
	func = Piecewise((0, True))
	def __init__(self, u, func = None):
		Xmu.__init__(self, u)
		if func is not None:
			self.set_xequals(func)


class UpwardGradientXmu(Xmu):
	""" UpwardGradientXmu provides the basis for creating upward straight-line membership functions, and generates the X-mu function automatically. """
	a = 0.0
	b = 0.0
	
	def __init__(self, u, a = None, b = None):
		""" Initialises the UpwardGradientXmu, and prepares its superclass (Xmu).
		@param	a	is the point where mu is no longer 0.0
		@param	b	is the point where mu becomes 1.0
		@returns	an instance of UpwardGradientXmu
		"""
		Xmu.__init__(self, u)
		if a is not None and b is not None:
			self.setMuFunction(float(a), float(b))
			self.setXFunction(float(a), float(b))
		
	def setAB(self, a, b):
		""" Setter for A and B, and (re)generates mu function and X-mu function in sympy.
		@param	a	is the point where mu is no longer 0.0
		@param	b	is the point where mu becomes 1.0
		@returns	a tuple of the results from setMuFunction and setXFunction. This is for simplicity.
		"""
		return self.setMuFunction(float(a), float(b)), self.setXFunction(float(a), float(b))
	
	def setMuFunction(self, a, b):
		""" Setter for the mu function based on a and b points.
		@param	a	is the point where mu is no longer 0.0
		@param	b	is the point where mu becomes 1.0
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		"""
		self.a = a
		self.b = b
		p = Piecewise(
			(((self.x - self.a)/(self.b - self.a)), ((self.a < self.x) & (self.x < self.b) )),
			(1.0, (self.x >= self.b)),
			(0.0, True)
		)
		self.set_muequals(p)
		return p

	def setXFunction(self, a, b):
		""" Setter for the X-mu function based on a and b points.
		@param	a	is the point where mu is no longer 0.0
		@param	b	is the point where mu becomes 1.0
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		"""
		self.a = a
		self.b = b
		x = Symbol('x', real=True)
		i = Interval(((ALPHA*self.b) - (ALPHA*self.a) + self.a), self.u.sup)
		
		self.set_xequals(i)
		return i


class DownwardGradientXmu(Xmu):
	""" DownwardGradientXmu provides the basis for creating downward straight-line membership functions, and generates the X-mu function automatically. """
	
	a = 0.0
	b = 0.0
	
	def __init__(self, u, a = None, b = None):
		""" Initialises the DownwardGradientXmu, and prepares its superclass (Xmu).
		@param	a	is the last point where mu is 1.0
		@param	b	is the point where mu becomes 0.0
		@returns	an instance of DownwardGradientXmu
		"""
		Xmu.__init__(self, u)
		if a is not None and b is not None:
			self.setMuFunction(float(a), float(b))
			self.setXFunction(float(a), float(b))
		
	def setAB(self, a, b):
		""" Setter for A and B, and (re)generates mu function and X-mu function in sympy.
		@param	a	is the last point where mu is 1.0
		@param	b	is the point where mu becomes 0.0
		@returns	a tuple of the results from setMuFunction and setXFunction. This is for simplicity.
		"""
		return self.setMuFunction(float(a), float(b)), self.setXFunction(float(a), float(b))
	
	def setMuFunction(self, a, b):
		""" Setter for the mu function based on a and b points.
		@param	a	is the last point where mu is 1.0
		@param	b	is the point where mu becomes 0.0
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		"""
		self.a = a
		self.b = b
		p = Piecewise(
			(((self.b - self.x)/(self.b - self.a)), ((self.a < self.x) & (self.x < self.b) )),
			(1.0, (self.x <= self.a)),
			(0.0, True)
		)
		self.set_muequals(p)
		return p

	def setXFunction(self, a, b):
		""" Setter for the X-mu function based on a and b points.
		@param	a	is the last point where mu is 1.0
		@param	b	is the point where mu becomes 0.0
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		"""
		self.a = a
		self.b = b
		x = Symbol('x', real=True)

		f = ((-ALPHA*self.b)+(ALPHA*self.a)+self.b)
		i = Interval(self.u.inf, f)
		
		self.set_xequals(i)
		return i

class TrapezoidalXmu(Xmu):
	""" TrapezoidalXmu provides the basis for creating Trapezoidal straight-line membership functions, and generates the X-mu function automatically. """
	
	a = 0.0
	b = 0.0
	c = 0.0
	d = 0.0
	
	def __init__(self, u, a = None, b = None, c = None, d = None):
		""" Initialises the TrapezoidalXmu, and prepares its superclass (Xmu).
		@param	a	is the last point where mu is 0.0
		@param	b	is the point where mu becomes 1.0
		@param	c	is the last point where mu is 1.0
		@param	d	is the point where mu becomes 0.0 again
		@returns	an instance of TrapezoidalXmu
		"""
		Xmu.__init__(self, u)
		if a is not None and b is not None and c is not None and d is not None:
			self.setMuFunction(float(a), float(b), float(c), float(d))
			self.setXFunction(float(a), float(b), float(c), float(d))
		
	def setAB(self, a, b, c, d):
		""" Setter for a,b,c,d , and (re)generates mu function and X-mu function in sympy.
		@param	a	is the last point where mu is 0.0
		@param	b	is the point where mu becomes 1.0
		@param	c	is the last point where mu is 1.0
		@param	d	is the point where mu becomes 0.0 again
		@returns	a tuple of the results from setMuFunction and setXFunction. This is for simplicity.
		"""
		return self.setMuFunction(float(a), float(b), float(c), float(d)), self.setXFunction(float(a), float(b), float(c), float(d))
		

	def setMuFunction(self, a, b, c, d):
		""" Setter for the mu function based on a, b, c and d points.
		@param	a	is the last point where mu is 0.0
		@param	b	is the point where mu becomes 1.0
		@param	c	is the last point where mu is 1.0
		@param	d	is the point where mu becomes 0.0 again
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		"""
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		
		p = Piecewise(
			(((self.x - self.a)/(self.b - self.a)), ((self.a < self.x) & (self.x < self.b) )),
			(1.0, ((self.b <= self.x) & (self.x <= self.c))),
			(((self.d - self.x)/(self.d - self.c)), ((self.c < self.x) & (self.x < self.d) )),
			(0.0, True)
		)
		self.set_muequals(p)
		return p
		
		return None
	
	def setXFunction(self, a, b, c, d):
		""" Setter for the X-mu function based on a, b, c and d points.
		@param	a	is the last point where mu is 0.0
		@param	b	is the point where mu becomes 1.0
		@param	c	is the last point where mu is 1.0
		@param	d	is the point where mu becomes 0.0 again
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		"""
		self.a = a
		self.b = b
		self.c = c
		self.d = d
		
		x = Symbol('x', real=True)
		
		f = ((ALPHA*self.b) - (ALPHA*self.a) + self.a)
		g = ((-ALPHA*self.d) + (ALPHA*self.c) + self.d)
		
		i = Interval(f, g)
		
		self.set_xequals(i)
		return i

class TriangularXmu(Xmu):
	""" TriangularXmu provides the basis for creating TriangularXmu straight-line membership functions, and generates the X-mu function automatically. """
	
	a = 0.0
	b = 0.0
	c = 0.0
	
	def __init__(self, u, a = None, b = None, c = None):
		""" Initialises the TriangularXmu, and prepares its superclass (Xmu).
		@param	a	is the last point where mu is 0.0
		@param	b	is the point where mu becomes 1.0
		@param	c	is the point where mu becomes 0.0 again
		@returns	an instance of TriangularXmu
		"""
		Xmu.__init__(self, u)
		if a is not None and b is not None and c is not None:
			self.setMuFunction(float(a), float(b), float(c))
			self.setXFunction(float(a), float(b), float(c))
	
	def setMuFunction(self, a, b, c):
		""" Setter for the mu function based on a, b and c points.
		@param	a	is the last point where mu is 0.0
		@param	b	is the point where mu becomes 1.0
		@param	c	is the point where mu becomes 0.0 again
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		@todo	clean-up, as this is currently a direct c&p of the Trapezoidal version
		"""
		self.a = a
		self.b = b
		self.c = b
		self.d = c
		
		p = Piecewise(
			(((self.x - self.a)/(self.b - self.a)), ((self.a < self.x) & (self.x < self.b) )),
			(1.0, ((self.b <= self.x) & (self.x <= self.c))),
			(((self.d - self.x)/(self.d - self.c)), ((self.c < self.x) & (self.x < self.d) )),
			(0.0, True)
		)
		self.set_muequals(p)
		return p
	
	def setXFunction(self, a, b, c):
		""" Setter for the X-mu function based on a, b and c points.
		@param	a	is the last point where mu is 0.0
		@param	b	is the point where mu becomes 1.0
		@param	c	is the point where mu becomes 0.0 again
		@returns	a sympy interval. This is for simplicity (so we don't have to call get.
		"""
		self.a = a
		self.b = b
		self.c = c
		
		x = Symbol('x', real=True)
		
		f = ((ALPHA*self.b) - (ALPHA*self.a) + self.a)
		g = ((-ALPHA*self.c) + (ALPHA*self.b) + self.c)
		
		i = Interval(f, g)
		
		self.set_xequals(i)
		return i

class Graph:
	""" Graph class extracts the plotting requirements away from the interface. So no need to know matplotlib, and could act as a unifying wrapper in the future.
	"""
	granularity = 100
	alphas = []
	plt = None
	u = Interval(0.0,1.0)
	
	def __init__(self, granularity, u):
		""" Initialises plotting library, and sets the granularity and universe.
		@param	granularity	the 'resolution' of the generated graph
		@param	u	the universe of the domain (a sympy interval)
		@return	an instantiated Graph object
		"""
		import matplotlib.pyplot as plt
		self.plt = plt
		self.set_granularity(granularity)
		self.set_u(u)
	
	def set_u(self, u):
		""" Sets u, the universe.
		@param	u	the universe
		"""
		self.u = u
	
	def set_plt(self, plt):
		""" Sets plt.
		@param	plt	in this version of the library, this is a matplotlib plt (plot) object.
		"""
		self.plt = plt
	
	def set_granularity(self, granularity):
		""" Sets the precision for graphing, and prepares generic alphas.
		@param	granularity
		"""
		self.granularity = granularity
		self.alphas = [float(i)/float(self.granularity) for i in range(0, granularity+1)]
		
	def prepare_plot(self, plot, plot_title, ylim_top=None):
		""" Sets up initial information for a plot.
		@param	plot	an Xmu instance to plot
		@param	plot_title	the title of the Xmu object
		@param	ylim_top	the maximum limit of the y axis
		"""
		self.add_plot(plot, plot_title)
		self.plt.xlabel('mu')
		self.plt.xlim(0.0, 1.0)
		if ylim_top == None:
			self.plt.ylim(float(self.u.inf), float(self.u.sup))
		else:
			self.plt.ylim(float(self.u.inf), float(ylim_top))
		self.plt.ylabel('X')
		
	
	def add_plot(self, plot, plot_title, colour=None):
		""" adds an Xmu object to the plot.
		@param	plot	an Xmu object
		@param	plot_title	the title of the Xmu object
		@param	colour	the colour of the visualisation of this Xmu object
		"""
		if isinstance(plot, Xmu):
			plot = plot.get_xequals()
		
		if colour is None:
			colour = rand(3,1)
		
		plot_xss_length = []
		plot_xss_start = []
		plot_alphas = []
		for i in self.alphas:
			subs = plot.subs(ALPHA, i)
			if not subs.is_EmptySet:
				if type(subs.args[0]) is not Interval:
					plot_xss_length.append(subs.sup - subs.inf)
					plot_xss_start.append(subs.inf)
					plot_alphas.append(i)
				else: # type(subs.args[0]) is Interval
					for part in subs.args:
						plot_xss_length.append(part.sup - part.inf)
						plot_xss_start.append(part.inf)
						plot_alphas.append(i)
		
		if(len(plot_alphas) > 0):
			self.plt.bar(left=plot_alphas, height=plot_xss_length, bottom=plot_xss_start, width=1.0/self.granularity, alpha=0.3, linewidth=0, label=plot_title, color=colour)
		else:
			print "Nothing to plot"
	
	
	
	def prepare_arithmetic_plot(self, A_interval, B_interval, func, plot_title, ylim_top=None):
		""" Prepares a purely arithmetic plot.
		@param	A_interval	one interval to add to a plot.
		@param	B_interval	another interval to add to a plot.
		@param	func	the arithmetic operation to apply on A using B.
		@param	plot_title	the title of this operation.
		@param	ylim_top	the top limit for the y axis.
		"""
		self.add_arithmetic_plot(A_interval, B_interval, func, plot_title)
		self.plt.xlabel('mu')
		self.plt.xlim(0.0, 1.0)
		if ylim_top == None:
			self.plt.ylim(float(self.u.inf), float(self.u.sup))
		else:
			self.plt.ylim(float(self.u.inf), float(ylim_top))
		self.plt.ylabel('X')
	
	
	
	def add_arithmetic_plot(self, A_interval, B_interval, func, plot_title, colour=None):
		""" Adds an arithmetic operation to the plot.
		@param	A_interval	one interval to add to a plot.
		@param	B_interval	another interval to add to a plot.
		@param	func	the arithmetic operation to apply on A using B.
		@param	plot_title	the title of this operation.
		@param	colour	the colour of this plot.
		"""
		if colour is None:
			colour = rand(3,1)
		
		plot_xss_length = []
		plot_xss_start = []
		plot_alphas = []
		for i in self.alphas:
			result = eval("A_interval." + func + "(B_interval, i)")
			if result is not None and not result.get_xequals().is_EmptySet:
				result = result.get_xequals()
				if type(result.args[0]) is not Interval:
					plot_xss_length.append(float(result.sup) - float(result.inf))
					plot_xss_start.append(float(result.inf))
					plot_alphas.append(float(i))
				else: # type(subs.args[0]) is Interval
					for part in result.args:
						plot_xss_length.append(float(part.sup) - float(part.inf))
						plot_xss_start.append(float(part.inf))
						plot_alphas.append(float(i))
		
		self.plt.bar(left=plot_alphas, height=plot_xss_length, bottom=plot_xss_start, width=1.0/float(self.granularity), alpha=0.3, linewidth=0, label=plot_title, color=colour)
	
	def show_plot(self):
		""" Shows a plot to screen. """
		self.plt.legend(loc='best')
		self.plt.show()



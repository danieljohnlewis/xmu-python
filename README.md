# The X-mu Python Library

## About
The X-mu Python Library implements the X-mu method of Fuzzy Set Theory using the Python Programming Language.

### Briefly about Fuzzy and X-mu
Fuzzy Set Theory is a little like classical set theory, where there are items of a set (or list). However, where it differs is that Fuzzy Sets allow items to be partially in/out of a set. This level of partiality is called the membership value. A membership value can be set either individually, with some kind of look-up table or basic formula, or through some kind of analytic membership function.

The X-mu methodology tries to solve some of the problems regarding traditional forms of fuzzy set theory, and it does this by inverting a membership function. The result allow easy analytic/symbolic merging of functions when performing set operations (such as intersection, union and difference) or fuzzy arithmetic operations (such as addition, subtraction, multiplication, division or power). This is where The X-mu Python Library comes in, it implements this methodology in the python programming language, and uses the power of SymPy (the symbolic computation engine for Python) to do this!

### Requirements
In order to run the X-mu Python Library the following is required:
* [Python 2](http://www.python.org/) - tested with version 2.7.6
* [SymPy](http://sympy.org/) - tested with version 0.7.4
* [NumPy](http://www.numpy.org/) - tested with version 1.7.1
* [Matplotlib](http://matplotlib.org/) - tested with version 1.3.1

## Authorship
The X-mu Python Library was created by Daniel Lewis while researching towards a PhD in Engineering Mathematics at the Intelligent Systems Laboratory of the University of Bristol. Daniel began X-mu research in 2012, and created this library in 2013. It was finally open-sourced in 2014.

Persistant information can be found about Daniel at his website:
http://vanirsystems.com/
He can be reached by email at:
daniel@vanirsystems.com

The X-mu methodology was primarily created through the work of Prof. Trevor Martin (University of Bristol), and Daniel Lewis (University of Bristol). It is being to various applications (including ratings, traffic data) and as a boosting method to other techniques (including Fuzzy Association Rule Mining, and Fuzzy Concept Lattices). A big thank you goes to Prof. Trevor Martin.

## License
Copyright 2014 Daniel Lewis

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

See the License for the specific language governing permissions and
limitations under the License.

## Basic Example
Here follows a basic example of use:

	from sympy import * # For Symbolic Computation
	from xmu import * # Our X-mu Library (local)
	init_printing() # Turns on pretty printing
	
	x = Symbol('x', real=True, bounded=True) # Creates a real value
	u = Interval(1.0, 6.0) # Creates u (our Universe)
	
	### Large
	large_a = 3.0
	large_b = 5.0
	large = UpwardGradientXmu(u, large_a, large_b) # An example of an Upward Gradient (i.e. an open Fuzzy Membership Function)
	
	vlarge_a = 5.0
	vlarge_b = 6.0
	vlarge = UpwardGradientXmu(u, vlarge_a, vlarge_b) # Another upward gradient
	
	### Small
	small_a = 2.0
	small_b = 4.0
	small = DownwardGradientXmu(u, small_a, small_b) # Downward gradient
	
	### Medium: Two options, either Triangular or Trapezoidal
	medium_a = 1
	medium_b = 3
	medium_c = 4
	medium_d = 6
	medium = TrapezoidalXmu(u, medium_a, medium_b, medium_c, medium_d)
	#medium = TriangularXmu(u, medium_a, medium_b, medium_d)
	
	### This will perform: (Small union Large) - Medium
	testplot = small.unionX(large).differenceX(medium).get_xequals()
	
	### Graph showing testplot and large
	graph = Graph(100, u)
	graph.prepare_plot(testplot, u"(Small union Large) - Medium")
	graph.add_plot(large, u"Large")
	graph.show_plot()


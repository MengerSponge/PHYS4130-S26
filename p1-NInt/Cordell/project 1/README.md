Code.py is a python script that demonstrates several numeric integration techniques.  % $`sin(sqrt(100x))^{2}`$ 0 to 2.  explanation of theory

The first technique it demonstrates is the trapezoid rule. The trapezoid rule takes the average of the leftpoint and rightpoint rules. These methods incorporate the Reimann sum technique which splits the area under the curve into N rectangles and adds up the area of each rectangle. The height of each rectangle is typically taken at the left, right, or middle of the rectangle--- hence the names of the rules. As one might imagine, the estimate becomes more accurate at higher values of N. Specifically, as N increases the error of leftpoint or rightpoint decreases proportionally to 1/N and the error of midpoint or trapezoid decreases proportionally to 1/N^2.

Simpson's method is a weighted average of the midpoint and trapezoid rule where midpoint is weighted by 2/3 and the trapezoid is weighted by 1/3. It is the most accurate out of the methods so far. As N increases, the error reduces by a factor of 1/N^4.

The Guassian quadrature method first parameterizes an integral to go from -1 to 1.

To run the code.py script open the command prompt from the directory you have the it in:
1) conda activate [the name of your environment]
2) python code.py

% document implementation of algortihms
Here is a breakdown and explanations of important code segments:

```python
def leftpoint(f,a,b,N):
    mysum = 0
    h = f
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    A_array = h(x_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[N-1]      
    return mysum

def rightpoint(f,a,b,N):
    mysum = 0
    h = f
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    A_array = h(x_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[0]     
    return mysum 

def trapezoid(f, a, b, N):
    return 0.5* (leftpoint(f, a, b, N) + rightpoint(f, a, b, N))

def midpoint(f,a,b,N):        
    h = f
    mysum = 0
    x_array = np.linspace(a, b, N)
    w = (b-a)/(N-1)
    mid_array = [x + (w/2) for x in x_array]
    mid_array = np.array(mid_array)
    A_array = h(mid_array)*w
    mysum = np.sum(A_array)
    mysum = mysum - A_array[N-1] 
    return mysum   

def Simpson(f, a, b, N):
    return (1/3 * trapezoid(f, a, b, N)) + (2/3 * midpoint(f, a, b, N))

def quad(a, b, N):
    roots, weights = sp.special.roots_legendre(N)
    x = ((b-a)*roots/2)+(a+b)/2
    dx_over_du = 2/(b-a)
    return dx_over_du* np.sum(weights*sin(x))
```

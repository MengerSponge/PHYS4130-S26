# Code Log
This is for me to commit what I've done code-wise, since I don't want to risk uploading MCNP stuff to github.


I created this starting from 4/28, but I did some coding before then starting Thursday 4/23. Before that I replicated the Godiva sphere and played around with some other geometries.

### 4/28
I created a geometry of concentric cylindrical surfaces with alternating layers of Uranium and water. There seems to be an issue with what I'm calling the inside and/or outside of each surface, so I would like to figure that out before I modify the density of Uranium such that the total Uranium mass in the geometry is 18 kg and estimate $k_{eff}$.

### 4/29
I fixed the bug. The issue is I didn't define a cell inside the inner-most cylinder so it defaulted to void, and my neutrons were spawning at the origin causing them to die instantly as I have it set to not track neutrons in void.

After I fixed the bug I estimated $k_{eff}$ = 1.73 Then, I changed from

CZ: Infinitely long cylindrical surfaces centered on z-axis
--> RCC: "Right circular cylinders" which are a little more involved to define but I can limit their height.

Using the following geometry:
```
c CELL CARDS
10   100  -18.74  +1 -2      imp:n=1                    
20   200  -1.0    +2 -3      imp:n=1                    
30   100  -18.74  +3 -4      imp:n=1                    
40   0            +4         imp:n=0                    
50   200  -1.0    -1         imp:n=1

c SURFACE CARDS
1   RCC     0 0 -50     0 0 100     7.000                                
2   RCC     0 0 -50     0 0 100     14.000
3   RCC     0 0 -50     0 0 100     21.000
4   RCC     0 0 -50     0 0 100     28.000
```

I estimated $k_{eff}$ = 1.70

However, this uses 2885 kg of Uranium. So now I need to shrink the volume of Uranium cells + density of Uranium as to only use 18 kg.

First, I tried taking the density down from 18.74 g/$cm^{3}$ to 0.12 (I did the math and it gave me 18 kg total mass). This gave me $k_{eff}$ = 0.57.

Then, I took the radius down some and bumped up the density. After doing a few rounds of adjustments plus simulations, here are some notable results:

For 9.74 kg total: This geometry I made yields $k_{eff}$ = 0.965. This is around the safety threshold and it only uses a little over half the Uranium they have in stock.

```
Practice 
c CELL CARDS
10   100  -18.74  +1 -2      imp:n=1                   
20   200  -1.0    +2 -3      imp:n=1                    
30   100  -18.74  +3 -4      imp:n=1                    
40   0            +4         imp:n=0                   
50   200  -1.0    -1         imp:n=1

c SURFACE CARDS
1   RCC     0 0 -50     0 0 100     10.000                                   
2   RCC     0 0 -50     0 0 100     10.500
3   RCC     0 0 -50     0 0 100     20.500
4   RCC     0 0 -50     0 0 100     21.000

c DATA CARDS
kcode 10000  1.0  100  200                    
ksrc  0.0  0.0  0.0                           
m100  92235 -.9473                             
      92238 -.0527
m200  1001   2 
      8016   1
```

For 16.1 kg total: This geometry I made yields $k_{eff}$ = 1.06, which achieves the goal of going critical with 18 kg of Uranium.

```
c CELL CARDS
10   100  -18.74  +1 -2      imp:n=1                 
20   200  -1.0    +2 -3      imp:n=1            
30   100  -18.74  +3 -4      imp:n=1        
40   0            +4         imp:n=0       
50   200  -1.0    -1         imp:n=1

c SURFACE CARDS
1   RCC     0 0 -40     0 0 80     10.000                                  
2   RCC     0 0 -40     0 0 80     11.000
3   RCC     0 0 -40     0 0 80     21.000
4   RCC     0 0 -40     0 0 80     22.000

c DATA CARDS
kcode 10000  1.0  100  200                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.9473               
      92238 -.0527
m200  1001   2 
      8016   1
```

### 4/30

See what happens when I add more layers of decreasing thickness with the ith layer from the center.

I changed the Uranium mix to be low-enriched 80/20
```
m100  92235 -.8000                              
      92238 -.2000
```
Now, $k_{eff}$ is 1.03 rather than 1.06.

I've been messing around with adding tops to the cylinders or leaving them open. I've found it seems more ideal to leave them open, so I can make them longer without having to worry about using more mass. 

I have been working with more layers but haven't been able to get a good $k_{eff}$ with more layers yet.

I realized that unfortunately I had previously not reached $k_{eff}$ with 18 kg because I forgot to account for the density 18.74 g/$cm^{3}$. So, I was actually using 16.1 x 18.74 kg = 302 kg. I have now accounted for the density correctly and I'm struggling to achieve criticality. My closest attempt so far is this geometry:

```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   0            +6         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -5.000     0 0 10.0     10.000
2   RCC     0 0 -5.050     0 0 10.1     10.100
3   RCC     0 0 -10.05     0 0 20.1     20.100
4   RCC     0 0 -10.10     0 0 20.2     20.200
5   RCC     0 0 -15.10     0 0 30.2     30.200
6   RCC     0 0 -15.15     0 0 30.3     30.300

c DATA CARDS
kcode 10000  1.0  100  200                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
which yields $k_{eff}$ = 0.72


Then,
```
1   RCC     0 0 -10     0 0 20.0     6.000                                  
2   RCC     0 0 -10     0 0 20.0     6.100
3   RCC     0 0 -10     0 0 20.0     12.100
4   RCC     0 0 -10     0 0 20.0     12.200
5   RCC     0 0 -10     0 0 20.0     15.200
6   RCC     0 0 -10     0 0 20.0     15.300
7   RCC     0 0 -10     0 0 20.0     18.300
8   RCC     0 0 -10     0 0 20.0     18.400
```
yields 0.65 with about 12 kg.

Lastly, 
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   0            +8         imp:n=0       

c SURFACE CARDS
c 1   RCC     0 0 -5.000     0 0 10.0     10.000                                  
c 2   RCC     0 0 -5.050     0 0 10.1     10.100
c 3   RCC     0 0 -10.05     0 0 20.1     20.100
c 4   RCC     0 0 -10.10     0 0 20.2     20.200
c 5   RCC     0 0 -15.10     0 0 30.2     30.200
c 6   RCC     0 0 -15.15     0 0 30.3     30.300
1   RCC     0 0 -20     0 0 40.0     5.000                                  
2   RCC     0 0 -20     0 0 40.0     5.100
3   RCC     0 0 -20     0 0 40.0     9.100
4   RCC     0 0 -20     0 0 40.0     9.200
5   RCC     0 0 -20     0 0 40.0     12.200
6   RCC     0 0 -20     0 0 40.0     12.300
7   RCC     0 0 -20     0 0 40.0     14.300
8   RCC     0 0 -20     0 0 40.0     14.400

c DATA CARDS
kcode 10000  1.0  100  200                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
yields $k_{eff}$ = 0.82 with about 18 kg.


Then,
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   0            +8         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -19     0 0 38.0     5.000                                  
2   RCC     0 0 -19     0 0 38.0     5.100
3   RCC     0 0 -19     0 0 38.0     9.100
4   RCC     0 0 -19     0 0 38.0     9.200
5   RCC     0 0 -19     0 0 38.0     12.200
6   RCC     0 0 -19     0 0 38.0     12.300
7   RCC     0 0 -19     0 0 38.0     15.300
8   RCC     0 0 -19     0 0 38.0     15.400

c DATA CARDS
kcode 10000  1.0  100  200                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
yields $k_{eff}$ = 0.85 with slightly less than 18 kg.

Then,
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   0            +8         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -18.5     0 0 37.0     4.000                                  
2   RCC     0 0 -18.5     0 0 37.0     4.100
3   RCC     0 0 -18.5     0 0 37.0     8.100
4   RCC     0 0 -18.5     0 0 37.0     8.200
5   RCC     0 0 -18.5     0 0 37.0     12.200
6   RCC     0 0 -18.5     0 0 37.0     12.300
7   RCC     0 0 -18.5     0 0 37.0     16.300
8   RCC     0 0 -18.5     0 0 37.0     16.400

c DATA CARDS
kcode 10000  1.0  100  200                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
yields $k_{eff}$ = 0.88 with slightly less than 18 kg

Then, 
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   200  -1.0    +8 -9      imp:n=1
100  100  -18.74  +9 -10      imp:n=1
110   0           +10         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -18     0 0 36.0     4.000                                  
2   RCC     0 0 -18     0 0 36.0     4.070
3   RCC     0 0 -18     0 0 36.0     8.070
4   RCC     0 0 -18     0 0 36.0     8.140
5   RCC     0 0 -18     0 0 36.0     12.140
6   RCC     0 0 -18     0 0 36.0     12.210
7   RCC     0 0 -18     0 0 36.0     16.210
8   RCC     0 0 -18     0 0 36.0     16.280
9   RCC     0 0 -18     0 0 36.0     20.280
10  RCC     0 0 -18     0 0 36.0     20.350

c DATA CARDS
kcode 10000  1.0  100  200                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
yields $k_{eff}$ = 0.92 for slightly less than 18kg.

This geomeotry,
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   200  -1.0    +8 -9      imp:n=1
100   0           +9         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -21     0 0 42.0     4.000                                  
2   RCC     0 0 -21     0 0 42.0     4.060
3   RCC     0 0 -21     0 0 42.0     8.060
4   RCC     0 0 -21     0 0 42.0     8.120
5   RCC     0 0 -21     0 0 42.0     12.120
6   RCC     0 0 -21     0 0 42.0     12.180
7   RCC     0 0 -21     0 0 42.0     16.180
8   RCC     0 0 -21     0 0 42.0     16.240
9   RCC     0 0 -21     0 0 42.0     50.000

c DATA CARDS
kcode 10000  1.0  100  200                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
yields $k_{eff}$ = 0.96 with about 12 kg of LEU.

Then, this geometry
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   200  -1.0    +8 -9      imp:n=1
100   0           +9         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -31     0 0 62.0     4.000                                  
2   RCC     0 0 -31     0 0 62.0     4.060
3   RCC     0 0 -31     0 0 62.0     8.060
4   RCC     0 0 -31     0 0 62.0     8.120
5   RCC     0 0 -31     0 0 62.0     12.120
6   RCC     0 0 -31     0 0 62.0     12.180
7   RCC     0 0 -31     0 0 62.0     16.180
8   RCC     0 0 -31     0 0 62.0     16.240
9   RCC     0 0 -31     0 0 62.0     50.000

c DATA CARDS
kcode 10000  1.0  20  50                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
yields $k_{eff}$ = 1.016251  std dev = 0.001535 with 18 kg LEU.

Then, this geometry:
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   200  -1.0    +8 -9      imp:n=1
100   0           +9         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -23     0 0 56.0     4.000                                  
2   RCC     0 0 -23     0 0 56.0     4.065
3   RCC     0 0 -23     0 0 56.0     8.065
4   RCC     0 0 -23     0 0 56.0     8.130
5   RCC     0 0 -23     0 0 56.0     12.130
6   RCC     0 0 -23     0 0 56.0     12.195
7   RCC     0 0 -23     0 0 56.0     16.195
8   RCC     0 0 -23     0 0 56.0     16.260
9   RCC     0 0 -23     0 0 56.0     50.000

c DATA CARDS
kcode 10000  1.0  20  50                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
yields $k_{eff}$ = 1.021767, std dev = 0.001480 with slightly less than 18 kg.

Then, adding more water on top:
```
c CELL CARDS
10   200  -1.0    -1         imp:n=1
20   100  -18.74  +1 -2      imp:n=1                 
30   200  -1.0    +2 -3      imp:n=1            
40   100  -18.74  +3 -4      imp:n=1     
50   200  -1.0    +4 -5      imp:n=1
60   100  -18.74  +5 -6      imp:n=1
70   200  -1.0    +6 -7      imp:n=1
80   100  -18.74  +7 -8      imp:n=1
90   200  -1.0    +8 -9      imp:n=1
100   0           +9         imp:n=0       

c SURFACE CARDS
1   RCC     0 0 -23     0 0 56.0     4.000                                  
2   RCC     0 0 -23     0 0 56.0     4.065
3   RCC     0 0 -23     0 0 56.0     8.065
4   RCC     0 0 -23     0 0 56.0     8.130
5   RCC     0 0 -23     0 0 56.0     12.130
6   RCC     0 0 -23     0 0 56.0     12.195
7   RCC     0 0 -23     0 0 56.0     16.195
8   RCC     0 0 -23     0 0 56.0     16.260
9   RCC     0 0 -50     0 0 100.0     50.000

c DATA CARDS
kcode 10000  1.0  20  50                     
ksrc  0.0  0.0  0.0                       
m100  92235 -.20               
      92238 -.80
m200  1001   2 
      8016   1
```
 final k(col/abs/trk len) = 1.037970     std dev = 0.000915
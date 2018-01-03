MCNP6 leakage study radius:10 cm, height:50 cm. 
c Cell Card
1 1 -6.25 -1 -2 3 imp:n=1  $ cylindrical reactor
2 0 2:-3:(1 -2 3) imp:n=0  $ outside world
 
c Surface Card
1 CZ 10
2 PZ 25.0
3 PZ -25.0

c Data Card
m1
     6012 -1.7806e-03
     6013 -2.0869e-05
     7014 -4.1485e-02
     8016 -3.5986e-03
     8017 -1.4569e-06
     74180 -2.4821e-07
     74182 -5.5423e-05
     74183 -3.0093e-05
     74184 -6.4787e-05
     74186 -6.0769e-05
     92235 -8.3666e-01
     92238 -1.1624e-01
kcode 150000 1 35 150
ksrc 0 0 0
     1 1 1
     -1 -1 -1 
     1 -1 1
     -1 1 1
     1 1 -1
mode n
print

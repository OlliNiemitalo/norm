Optimize polynomial coefficients for Muon Newton-Schulz iteration, using Differential Evolution.

This work is placed in the public domain / CC0.

Compile:

```shell
g++ optimize.cpp opti.cpp -ffast-math -march=native -O3 -Wno-unused-result
```

Run:

```shell
./a.out
```

## Version 1: no cushioning, no cumulative error

The code picks a desired number of x values for error (cost) evaluation, using not uniform distribution over the requested interval but by adding more points near the ends, similar to Chebyshev nodes. In the beginning of the optimization, mean square error is used as a proxy for max abs error and between desired steps the error used is ramped linearly to max abs error. The coefficient of the linear term is shared between all polynomials. Otherwise the results might differ in a redundant way between optimization runs.

The current optimization problem parameters are: 8192 sample points, ramp between steps 500k and 1M (ramp will take place well before optimization convergence), 5 polynomials, 5th degree polynomials, optimization x interval 0.001 to 1.010916328. The current Differential Evolution parameters are: cross-over 0.999, weight of difference 0.76, population size 1000.

Current best result for x = [0.001, 1.010916328], plotted in [Desmos](https://www.desmos.com/calculator):

![image](https://github.com/user-attachments/assets/0806a8c9-04e8-42a1-ac64-88910a0bc7b5)

```
(4.01528731974196961829, -11.64682000026122210556, 8.45593805664630160379),
(4.01528731974196961829, -12.99585762597172688970, 10.56972369767224506631),
(4.01528731974196961829, -14.15392102240266680724, 12.74917448599765812389),
(4.01528731974196961829, -13.88019768698666744910, 13.15739447087707780781),
(4.01528731974196961829, -9.01530701517418542323, 7.06338829287047698102),

4.01528731974196961829 x^1 + -11.64682000026122210556 x^3 + 8.45593805664630160379 x^5
4.01528731974196961829 x^1 + -12.99585762597172688970 x^3 + 10.56972369767224506631 x^5
4.01528731974196961829 x^1 + -14.15392102240266680724 x^3 + 12.74917448599765812389 x^5
4.01528731974196961829 x^1 + -13.88019768698666744910 x^3 + 13.15739447087707780781 x^5
4.01528731974196961829 x^1 + -9.01530701517418542323 x^3 + 7.06338829287047698102 x^5

Best cost 0.115707
```

This can be compared to cost 0.1398750080072837 of Polar Express (https://arxiv.org/abs/2505.16932), found by printing `l` and `u` after each iteration of their code and subtracting 1 from the last `u`. I got my x range endpoints by numerically solving the x for which the composite function is equal to `l` or `u`.

```
l u
0.008205137512087042 1.991794862487913
0.033368021342595415 1.9666319786574047
0.13048059786589022 1.8695194021341097
0.4259904513958469 1.574009548604153
0.8601249919927162 1.1398750080072837

8.205160414005569*x^1 + -22.90193498705604*x^3 + 16.460724910180307*x^5
4.066915619879589*x^1 + -2.861284534588477*x^3 + 0.5183804464778609*x^5
3.9134926112054607*x^1 + -2.824251876723087*x^3 + 0.5248485625148532*x^5
3.3060139701337707*x^1 + -2.430227567449687*x^3 + 0.4869515205509482*x^5
2.304016813944474*x^1 + -1.6427206546268964*x^3 + 0.4009100949022211*x^5
```

Naively looking it's not optimal in minimum max abs error sense, but they do include a safety margin at every iteration rather than just in the beginning as I do (I guess). Intuitively their "greedy is optimal" theory seems correct. Comparison of my f(x) and their g(x) composite functions:

![image](https://github.com/user-attachments/assets/7bb54ee5-ce6f-4a27-aef8-de2157553f65)

![image](https://github.com/user-attachments/assets/bec0d51e-dec8-4c6d-add7-16dcb3047596)

![image](https://github.com/user-attachments/assets/49feba92-a0de-4e0a-bf5c-e632215b2c6f)

Truncated compositions using my f(x) are equioscillating, with $2 \times 3^k - 2$ stationary points (counting also those for negative x) for each degree- $5^k$ composition of the $k$ first degree-5 polynomials:

![image](https://github.com/user-attachments/assets/af3061b2-43e9-4738-872f-a8925b06ac55)

## With cumulative error

After including a multiplicative error of 1.01 at every iteration, the maximum absolute deviation cost became a bit higher, 0.12710828643358684786. I also removed the RMS error and ramping trick, for simplicity, and still get convergence. I also enabled early-out in the cost function, which now terminates early if the solution proposal won't be used by the evolutionary algorithm based on a lower limit on the cost for the solution proposal.

```
(3.99971006509232740456, -11.85609455039494619655, 8.79678969844277425238),
(3.99971006509232740456, -13.06980619978764934785, 10.73203704365077726379),
(3.99971006509232740456, -14.06685173761908735912, 12.63933777733973329305),
(3.99971006509232740456, -13.66404203747230461374, 12.77831792129320476192),
(3.99971006509232740456, -8.89224179372500422858, 6.84907182425005967019),

3.99971006509232740456 x^1 + -11.85609455039494619655 x^3 + 8.79678969844277425238 x^5
3.99971006509232740456 x^1 + -13.06980619978764934785 x^3 + 10.73203704365077726379 x^5
3.99971006509232740456 x^1 + -14.06685173761908735912 x^3 + 12.63933777733973329305 x^5
3.99971006509232740456 x^1 + -13.66404203747230461374 x^3 + 12.77831792129320476192 x^5
3.99971006509232740456 x^1 + -8.89224179372500422858 x^3 + 6.84907182425005967019 x^5

Best cost 0.127108
```

I also tried to include something similar to their "cushion", but that led to promoting exceedingly large coefficients in the polynomials so I ditched that effort.
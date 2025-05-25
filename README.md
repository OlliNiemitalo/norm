Optimize polynomial coefficients for Muon Newton-Schulz iteration, using Differential Evolution.

Compile:

```shell
g++ optimize.cpp opti.cpp -ffast-math -march=native -O3 -Wno-unused-result
```

Run:

```shell
./a.out
```

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

Some coefficients are quite big! And the intermediate values go to almost zero. Numerically it's not good.

This can be compared to cost 0.1398750080072837 of Polar Express (https://arxiv.org/abs/2505.16932), found by printing `l` and `u` after each iteration of their code and subtracting 1 from the last `u`. I got my x range endpoints by numerically solving the x for which the composite function is `l` and `u`.

```
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

I don't know if their result was affected by their numerical safety considerations but it's not optimal in minimum max abs error sense.
